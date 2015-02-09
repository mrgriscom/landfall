from Queue import Queue
import proc
import geodesy
import math
from shapely.geometry import Polygon, MultiPolygon
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


def rectify_poly(bound, clockwise=True, rollover=180):
    # fix if crosses IDL or contains a pole. will not contain holes, though
    _ix = lambda ix: ix%len(bound)
    _next = lambda i: _ix(i+1)

    crossings = {}
    for i in xrange(len(bound)):
        p0 = bound[i]
        p1 = bound[_next(i)]
        if abs(p0[0] - p1[0]) > 180.:
            if p0[0] > p1[0]:
                left, right = p0, p1
            else:
                left, right = p1, p0
            to_idl = rollover - left[0]
            lon_diff = (right[0] - left[0]) % 360.
            crossing = left[1] + (right[1] - left[1]) * to_idl / lon_diff
            crossings[_next(i)] = crossing
    if not crossings:
        # this ignores the 'negative hole' scenario (dir of poly doesn't match 'clockwise')
        return Polygon(bound)

    def crossing_point(lat, ref):
        return (rollover - (360 if abs(rollover - ref[0]) > 180. else 0), lat)

    segments = {}
    for i, lat in crossings.iteritems():
        ix = i
        seg = [crossing_point(lat, bound[ix])]
        while True:
            seg.append(bound[ix])
            ix = _next(ix)
            if ix in crossings:
                break
        seg.append(crossing_point(crossings[ix], bound[_ix(ix-1)]))
        segments[i] = (seg, ix)

    crossings_by_lat = sorted(crossings.keys(), key=lambda i: crossings[i])
    ix_xbl = dict((e, i) for i, e in enumerate(crossings_by_lat))
    bounds = []
    while segments:
        start = None
        ix = None
        while True:
            if ix is None:
                ix, (seg, following_ix) = segments.popitem()
                start = ix
                subbound = list(seg)
            else:
                seg, following_ix = segments.pop(ix)
                subbound.extend(seg)

            end_right = (seg[-1][0] == rollover)
            dir_up = (end_right != clockwise)

            next_ix_ix = ix_xbl[following_ix] + (1 if dir_up else -1)
            if next_ix_ix >= 0 and next_ix_ix < len(crossings_by_lat):
                next_ix = crossings_by_lat[next_ix_ix]
            else:
                pole_lat = 90 if dir_up else -90
                pole_patch = [(rollover - 360., pole_lat), (rollover, pole_lat)]
                if end_right:
                    pole_patch.reverse()
                subbound.extend(pole_patch)
                next_ix = following_ix

            if next_ix == start:
                break
            ix = next_ix
        bounds.append(subbound)

    return MultiPolygon([Polygon(b) for b in bounds])

def swap(p):
    return (p[1], p[0])
def antipode(p):
    return (-p[0], p[1]+180.)

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
        seg = [swap(p) for p in project_line(s, 10.)[:-1]]
        if bound and geodesy.distance(swap(bound[-1]), swap(seg[-1])) < .001:
            continue
        bound.extend(seg)
    return rectify_poly(bound, bearing1 > bearing0)

def project_line(func, tolerance, maxlen=60.):
    line = []

    def _proj(a, b, fa, fb):
        lonb = geodesy.anglenorm(fb[1], 180-fa[1])
        fmid_interp = [
            .5*(fa[0] + fb[0]),
            geodesy.anglenorm(.5*(fa[1] + lonb)),
        ]
        mid = .5 * (a + b)
        fmid = func(mid)

        def error(tolerance):
            return geodesy.distance(fmid, fmid_interp) > tolerance
        def linlen(maxlen):
            lindist = ((fa[0]-fb[0])**2. + (fa[1]-lonb)**2.)**.5
            if lindist > maxlen:
                geodist = geodesy.distance(fa, fb)
                if geodist > 100.:
                    return True
            return False
        split = (error(tolerance) or linlen(maxlen))

        if split:
            #print a, b, fa, fb, fmid, fmid_interp, error(), linlen()
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
    # test if inside land before starting?
    minmerc = dist_to_merc(near_dist)
    maxmerc = dist_to_merc(.5*EARTH_CIRCUMF)
    zeromerc = dist_to_merc(0.)
    result = viewshed_drilldown(ix, p, bearing_res, bearing0, bearing1, minmerc, maxmerc, MaskNode())
    # and past antipode... (think bearing offset will break for N/S poles)
    viewshed_drilldown(ix, antipode(p), bearing_res, 180.-bearing0, 180.-bearing1, zeromerc, maxmerc, result)
    return result

EARTH_CIRCUMF = geodesy.EARTH_MEAN_RAD * 2. * math.pi
def dist_to_merc(dist):
    MIN_DIST = 1.
    dist = min(max(dist, MIN_DIST), .5*EARTH_CIRCUMF - MIN_DIST)
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

    def dump(self, min, max, res=None):
        mid = .5*(min + max)
        span = max - min
        val = None
        if self.null():
            val = -1
        elif self.val is not None:
            val = self.val

        if val is not None:
            if res is None or span <= res:
                yield (mid, span, val)
            else:
                for e in self.dump(min, mid, res):
                    yield e
                for e in self.dump(mid, max, res):
                    yield e
        else:
            for e in self.left.dump(min, mid, res):
                yield e
            for e in self.right.dump(mid, max, res):
                yield e

def viewshed_drilldown(ix, p, bearing_res, blo, bhi, dlo, dhi, mask):
    print blo, bhi, merc_to_dist(dlo), merc_to_dist(dhi)

    # bail if this region is already fully masked
    if mask.complete:
        print '  masked'
        return

    bspan = abs(bhi - blo)
    split_horiz, split_vert = False, False
    if bspan > 180.:
        split_horiz = True
    else:
        aspect = (dhi - dlo) / math.radians(bspan)
        if aspect > 2**.5 + 1e-6:
            split_vert = True
        elif aspect < 2**-.5 - 1e-6:
            split_horiz = True

    if split_horiz:
        print 'split-h'
        bmid = .5*(blo+bhi)
        mask.setLeft(viewshed_drilldown(ix, p, bearing_res, blo, bmid, dlo, dhi, mask.getLeft()))
        mask.setRight(viewshed_drilldown(ix, p, bearing_res, bmid, bhi, dlo, dhi, mask.getRight()))
        return mask
    if split_vert:
        print 'split-v'
        dmid = .5*(dlo+dhi)
        viewshed_drilldown(ix, p, bearing_res, blo, bhi, dlo, dmid, mask)
        viewshed_drilldown(ix, p, bearing_res, blo, bhi, dmid, dhi, mask)
        return mask

    view = make_view_poly(p, merc_to_dist(dlo), merc_to_dist(dhi), blo, bhi)
    hit = test_poly(ix, view)
    if hit == TYPE_WATER:
        print '  water'
        return

    if bspan <= bearing_res:
        print '  terminal %s' % hit
        dist = merc_to_dist(dlo)
        if bhi < blo:
            # past antipode
            dist += .5*EARTH_CIRCUMF
        mask.setVal(dist)
    else:
        print '  recurse'
        bmid = .5*(blo+bhi)
        dmid = .5*(dlo+dhi)
        mask.setLeft(viewshed_drilldown(ix, p, bearing_res, blo, bmid, dlo, dmid, mask.getLeft()))
        mask.setRight(viewshed_drilldown(ix, p, bearing_res, bmid, bhi, dlo, dmid, mask.getRight()))
        mask.setLeft(viewshed_drilldown(ix, p, bearing_res, blo, bmid, dmid, dhi, mask.getLeft()))
        mask.setRight(viewshed_drilldown(ix, p, bearing_res, bmid, bhi, dmid, dhi, mask.getRight()))
    return mask

def showpano(node, min, max):
    with open('/tmp/bbb', 'w') as f:
        for bear, width, dist in node.dump(min, max):
            if dist < 0:
                f.write('\n')
            else:
                f.write('%s %s\n' % (bear, -math.log(dist)))
    from subprocess import Popen, PIPE
    p = Popen('gnuplot', stdin=PIPE)
    p.stdin.write('plot "/tmp/bbb" using 1:2 with lines\n')

def tomap(node, ref, bmin, bmax):
    path = []
    prevdist = None
    antip = antipode(ref)
    for bear, width, dist in node.dump(bmin, bmax):
        p = geodesy.plot(ref, bear, dist if dist > 0 else (geodesy.EARTH_MEAN_RAD*math.pi - 1000))[0]
        if prevdist is not None and (dist > .5*EARTH_CIRCUMF) != (prevdist > .5*EARTH_CIRCUMF):
            path.append(antip)
        path.append(p)
        prevdist = dist

    def line(points):
        return {
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": [(p[1], p[0]) for p in points],
            },
            "properties": {
            }
        }
    SEGSIZE = 1000
    segments = [path[i:min(i+1+SEGSIZE, len(path))] for i in xrange(0, len(path), SEGSIZE)]
    segments.append([geodesy.plot(ref, bear, 10000)[0] for bear in xrange(0, 361, 5)])
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

def dump(node, p, minres, start, end, mindist):
    postings = list(node.dump(start, end, minres))
    res = postings[0][1]
    data = {
        'origin': p,
        'range': [start, end],
        'res': res,
        'min_dist': mindist,
        'postings': [k[2] for k in postings],
    }
    import tempfile
    path = tempfile.mktemp()
    with open(path, 'w') as f:
        json.dump(data, f, indent=True)
    return path

def _vs(ix, p, bearing_res, bearing0, bearing1, near_dist):
    result = viewshed(ix, p, bearing_res, bearing0, bearing1, near_dist)
    return {
        'result': result,
        'map': lambda: tomap(result, p, bearing0, bearing1),
        'plot': lambda: showpano(result, bearing0, bearing1),
        'dump': lambda: dump(result, p, bearing_res, bearing0, bearing1, near_dist),
    }
