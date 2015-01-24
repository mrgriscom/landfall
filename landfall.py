from Queue import Queue
import proc
import geodesy
import math
from shapely.geometry import Polygon
import json

OVL_NONE = 0
OVL_PARTIAL = 1
OVL_FULL = 2

TYPE_LAND = 1
TYPE_WATER = -1
TYPE_MIXED = 0

# worried about 10m resolution if 'near' boundary is close

def test_poly(ix, poly):
    # TODO trim poly on drilldown?
    # TODO return largest enclosing quad?

    fringe = Queue()
    terminal_types = set()

    fringe.put((ix, proc.ROOT))

    while not fringe.empty():
        node, ext = fringe.get()
        is_terminal = not hasattr(node, '__iter__')
        quad = proc.mkquad(ext)

        if quad.intersects(poly):
            if quad.within(poly):
                status = OVL_FULL
            else:
                status = OVL_PARTIAL
        else:
            status = OVL_NONE

        if status == OVL_NONE:
            continue
        elif is_terminal or status == OVL_FULL:
            if not is_terminal or node == TYPE_MIXED:
                return TYPE_MIXED
            terminal_types.add(node)
            if len(terminal_types) > 1:
                return TYPE_MIXED
        elif status == OVL_PARTIAL:
            for subnode, subext in zip(node, proc.quadchildren(ext)):
                fringe.put((subnode, subext))
    return terminal_types.pop()


def rectify_poly(poly):
    # fix if crosses IDL or contains a pole. will not contain holes, though
    return poly


def make_view_poly(p, near, far, wherever, you_are):
    bearing0, bearing1 = wherever, you_are
    assert (bearing0 % 360.) != (bearing1 % 360.), 'can\'t handle complete field of view'

    def greatcircle(bearing):
        pd = geodesy.plotter_d(p, bearing)
        return lambda t: pd((1-t) * near + t * far)[0]
    def distarc(dist):
        pb = geodesy.plotter_b(p, dist)
        return lambda t: pb((1-t) * bearing0 + t * bearing1)
    def rev(f):
        return lambda t: f(1-t)

    sides = [
        greatcircle(bearing0),
        distarc(far),
        rev(greatcircle(bearing1)),
        rev(distarc(near)),
    ]
    bound = []
    for s in sides:
        # TODO remove zero-length sides
        bound.extend((p[1], p[0]) for p in project_line(s, 10.)[:-1])
    return Polygon(bound)

def project_line(func, tolerance, maxlen=60.):
    line = []

    def _proj(a, b, fa, fb):
        fmid_interp = [
            .5*(fa[0] + fb[0]),
            geodesy.anglenorm(.5*(fa[1] + geodesy.anglenorm(fb[1], 180-fa[1]))), # handle crossing of IDL
        ]
        mid = .5 * (a + b)
        fmid = func(mid)

        error = lambda: geodesy.distance(fmid, fmid_interp)
        linlen = lambda: ((fa[0]-fb[0])**2. + (fa[1]-fb[1])**2.)**.5
        split = (error() > tolerance or linlen() > maxlen)

        if split:
            _proj(a, mid, fa, fmid)
            line.append(fmid)
            _proj(mid, b, fmid, fb)

    fa = func(0.)
    fb = func(1.)
    line.append(fa)
    _proj(0., 1., fa, fb)
    line.append(fb)
    return line

def merc(rad):
    return math.log(math.tan(math.pi/4 + rad/2))
def invmerc(mrad):
    return 2*math.atan(math.exp(mrad)) - math.pi/2

def viewshed(ix, p, bearing_res, bearing0, bearing1, near_dist):
    # TODO add dist resolution
    minmerc = dist_to_merc(near_dist)
    return viewshed_drilldown(ix, p, bearing_res, bearing0, bearing1, minmerc, dist_to_merc(5e6)+0*-minmerc, MaskNode())

def dist_to_merc(dist):
    return merc(dist/geodesy.EARTH_MEAN_RAD - math.pi/2)
def merc_to_dist(mrad):
    return (invmerc(mrad) + math.pi/2) * geodesy.EARTH_MEAN_RAD

class MaskNode(object):
    def __init__(self):
        self.val = None
        self.left = None
        self.right = None
        self.complete = False

    def null(self):
        return self.val is None and self.left is None and self.right is None

    def getLeft(self):
        if not self.left:
            self.left = MaskNode()
        return self.left
    def getRight(self):
        if not self.right:
            self.right = MaskNode()
        return self.right

    def setLeft(self, n):
        if n:
            self.left = n
            self.checkComplete()
    def setRight(self, n):
        if n:
            self.right = n
            self.checkComplete()
    def checkComplete(self):
        if all(n and n.complete for n in (self.left, self.right)):
            self.complete = True

    def setVal(self, val):
        if self.val is not None:
            return
        elif self.null():
            self.val = val
        else:
            self.left.setVal(val)
            self.right.setVal(val)
        self.complete = True

    def dump(self, min, max):
        mid = .5*(min + max)
        if self.null():
            yield (mid, -1)
        elif self.val is not None:
            yield (mid, self.val)
        else:
            for e in self.left.dump(min, mid):
                yield e
            for e in self.right.dump(mid, max):
                yield e

def viewshed_drilldown(ix, p, bearing_res, blo, bhi, dlo, dhi, mask):
    print blo, bhi, merc_to_dist(dlo), merc_to_dist(dhi)

    # bail if this region is already fully masked
    if mask.complete:
        print '  masked'
        return

    bspan = bhi - blo
    # split based on lonspan
    # split based on aspect ratio
    """
    if bearing_span > 180.:
        split_horiz
    aspect = (ymax - ymin) / math.radians(bearing_span)
    if aspect > 2**.5:
        split_vert
    elif aspect < 2**-.5:
        split horiz
    """

    view = make_view_poly(p, merc_to_dist(dlo), merc_to_dist(dhi), blo, bhi)
    hit = test_poly(ix, view)
    if hit == TYPE_WATER:
        print '  water'
        return

    if bspan <= bearing_res:
        print '  terminal %s' % hit
        mask.setVal(merc_to_dist(dlo))
    else:
        print '  recurse'
        bmid = .5*(blo+bhi)
        dmid = .5*(dlo+dhi)
        mask.setLeft(viewshed_drilldown(ix, p, bearing_res, blo, bmid, dlo, dmid, mask.getLeft()))
        mask.setRight(viewshed_drilldown(ix, p, bearing_res, bmid, bhi, dlo, dmid, mask.getRight()))
        mask.setLeft(viewshed_drilldown(ix, p, bearing_res, blo, bmid, dmid, dhi, mask.getLeft()))
        mask.setRight(viewshed_drilldown(ix, p, bearing_res, bmid, bhi, dmid, dhi, mask.getRight()))
    return mask


def render(ix, depth):
    dim = 2**depth
    BITMAP = [[0 for i in xrange(dim)] for j in xrange(dim)]
    def paint(val, z, (x, y)):
        c = {-1: 0, 0: 128, 1: 255}[val]
        celldim = 2**(depth - z)
        for i in xrange(celldim):
            for j in xrange(celldim):
                BITMAP[y*celldim + j][x*celldim + i] = c
    proc_index(ix, paint)

    import os
    import tempfile
    raw = tempfile.mktemp('.grey')
    img = tempfile.mktemp('.png')

    with open(raw, 'w') as f:
        f.write(''.join(''.join(chr(col) for col in row) for row in BITMAP))
    os.popen('convert -size %dx%d -depth 8 gray:%s %s' % (dim, dim, raw, img))
    os.popen('gnome-open %s' % img)




def poc(ix, p, near, far, b0, b1, x, y):
    BITMAP = [[0 for i in xrange(x)] for j in xrange(y)]

    deltad = (far - near) / float(y)
    deltab = (b1 - b0) / float(x)
    for j in xrange(y-1, -1, -1):
        for i in xrange(x):
            print i, j
            nr = near + deltad * j
            fr = near + deltad * (j + 1)
            left = b0 + deltab * i
            right = b0 + deltab * (i + 1)
            view = make_view_poly(p, nr, fr, left, right)

            res = test_poly(ix, view)
            BITMAP[y-1-j][i] = {-1: 0, 0: 128, 1: 255}[res]

    import os
    import tempfile
    raw = tempfile.mktemp('.grey')
    img = tempfile.mktemp('.png')
    with open(raw, 'w') as f:
        f.write(''.join(''.join(chr(col) for col in row) for row in BITMAP))
    os.popen('convert -size %dx%d -depth 8 gray:%s %s' % (x, y, raw, img))
    os.popen('gnome-open %s' % img)

# main algo



def showview(poly):
    with open('/tmp/aaa', 'w') as f:
        f.write('\n'.join('%s %s' % c for c in poly.exterior.coords))
    from subprocess import Popen, PIPE
    p = Popen('gnuplot', stdin=PIPE)
    p.stdin.write('plot "/tmp/aaa" using 2:1 with lines\n')

def showpano(node, min, max):
    with open('/tmp/bbb', 'w') as f:
        for bear, dist in node.dump(min, max):
            if dist < 0:
                f.write('\n')
            else:
                f.write('%s %s\n' % (bear, -math.log(dist)))
    from subprocess import Popen, PIPE
    p = Popen('gnuplot', stdin=PIPE)
    p.stdin.write('plot "/tmp/bbb" using 1:2 with lines\n')

def tomap(node, ref, bmin, bmax):
    path = []
    for bear, dist in node.dump(bmin, bmax):
        p = geodesy.plot(ref, bear, dist if dist > 0 else (geodesy.EARTH_MEAN_RAD*math.pi - 1000))[0]
        path.append((p[1], p[0]))

    def line(points):
        return {
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": points,
            },
            "properties": {
            }
        }
    SEGSIZE = 1000
    segments = [path[i:min(i+1+SEGSIZE, len(path))] for i in xrange(0, len(path), SEGSIZE)]
    data = {
        "type": "FeatureCollection",
        "features": map(line, segments),
    }
    with open('/tmp/ccc', 'w') as f:
        json.dump(data, f)
    import os
    os.popen('ogr2ogr -f kml /tmp/test.kml /tmp/ccc')
    with open('/tmp/test.kml') as f:
        kml = f.read()
    kml = kml.replace('<LineString><coordinates>', '<LineString><altitudeMode>clampToGround</altitudeMode><tessellate>1</tessellate><coordinates>')
    with open('/tmp/test.kml', 'w') as f:
        f.write(kml)
    os.popen('gnome-open /tmp/test.kml')
