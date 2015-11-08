from landfall import *
from shapely.geometry import LineString, LinearRing, MultiLineString, GeometryCollection

OVERLAP_NONE = 0
OVERLAP_PARTIAL = 1
OVERLAP_FULL = 2

MAX_DEPTH = 13  # ~5km, less than half of typical nautical claim (12 nmi)

class Quad(object):
    def __init__(self, x0, x1, y0, y1):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1

    @property
    def poly(self):
        return box(self.x0, self.y0, self.x1, self.y1)
        
    @staticmethod
    def root():
        return Quad(-180, 180, 90, -270)

    def children(self):
        xmid = .5 * (self.x0 + self.x1)
        ymid = .5 * (self.y0 + self.y1)
        yield Quad(self.x0, xmid, self.y0, ymid)
        yield Quad(xmid, self.x1, self.y0, ymid)
        yield Quad(self.x0, xmid, ymid, self.y1)
        yield Quad(xmid, self.x1, ymid, self.y1)

    def clip(self, polys):
        quad = self.poly
        subpolys = []
        for poly in polys:
            if quad.intersects(poly):
                if quad.within(poly):
                    return OVERLAP_FULL, None
                else:
                    if not quad.contains(poly):
                        subpoly = quad.intersection(poly)
                    else:
                        subpoly = poly
                    subpolys.append(subpoly)
        if subpolys:
            return OVERLAP_PARTIAL, subpolys
        else:
            return OVERLAP_NONE, None

def spatial_decompose(geometry, max_depth, verbose=False):
    ovl_full = set()
    ovl_partial = {}

    def descend(ix, quad, geom):
        if verbose:
            print (ix or '~'), len(geom)
        depth = len(ix)
        ovl, subgeom = quad.clip(geom)
        if ovl == OVERLAP_NONE:
            pass
        elif ovl == OVERLAP_FULL:
            ovl_full.add(ix)
        else:
            if depth == max_depth:
                ovl_partial[ix] = subgeom
            else:
                for i, subquad in enumerate(quad.children()):
                    descend(ix + str(i), subquad, subgeom)

    descend('', Quad.root(), geometry)
    return ovl_full, ovl_partial

def all_rings(polys):
    for poly in polys:
        for p in explode_poly(poly):
            yield p.exterior
            for r in p.interiors:
                yield r

def explode_lines(geom):
    if any(isinstance(geom, t) for t in (LineString, LinearRing)):
        return [geom]
    elif isinstance(geom, MultiLineString):
        return geom
    elif isinstance(geom, GeometryCollection):
        return [e for e in geom if any(isinstance(e, t) for t in (LineString, LinearRing))]
    else:
        return []

def render(decomp, depth):
    full, partial = decomp

    dim = 2**depth
    BITMAP = [[0 for i in xrange(dim)] for j in xrange(dim)]
    def paint(ix, val):
        z, x, y = 0, 0, 0
        for k in ix:
            k = int(k)
            z += 1
            x = 2*x + k % 2
            y = 2*y + k // 2
        c = {0: 0, 1: 128, 2: 255}[val]
        celldim = 2**(depth - z)
        for i in xrange(celldim):
            for j in xrange(celldim):
                BITMAP[y*celldim + j][x*celldim + i] = c
    for ix in full:
        paint(ix, OVERLAP_FULL)
    for ix in partial.keys():
        paint(ix, OVERLAP_PARTIAL)

    import os
    import tempfile
    raw = tempfile.mktemp('.grey')
    img = tempfile.mktemp('.png')

    with open(raw, 'w') as f:
        f.write(''.join(''.join(chr(col) for col in row) for row in BITMAP))
    os.popen('convert -size %dx%d -depth 8 gray:%s %s' % (dim, dim, raw, img))
    os.popen('gnome-open %s' % img)

def build_admin_index(path):
    admins = load_geometry(path, 'CODE')
    print '%d areas' % len(admins)

    ixfull = collections.defaultdict(set)
    ixpartial = collections.defaultdict(dict)

    for k, v in admins.items():
        full, partial = spatial_decompose(v, MAX_DEPTH)
        print k, len(full), len(partial)
        #render((full, partial), MAX_DEPTH)
        for e in full:
            ixfull[e].add(k)
        for e, f in partial.iteritems():
            ixpartial[e][k] = f

    return ixfull, ixpartial

def load_coastline(path):
    coastpolys = load_geometry(path)
    coastline = list(all_rings(coastpolys))
    full, partial = spatial_decompose(coastline, MAX_DEPTH, verbose=True)
    assert not full
    return partial

def ix_ancestor(ix):
    for i in xrange(len(ix), -1, -1):
        yield ix[:i]

def analyze(coast, admin):
    admin_full, admin_partial = admin
    partcount = collections.defaultdict(lambda: 0)
    fullonlycount = 0
    nullius = set()
    for ix, _ in coast.iteritems():
        if ix[0] in '23':
            continue

        ap = set(admin_partial.get(ix, {}).keys())
        af = reduce(lambda a, b: a.union(b), (admin_full.get(anc_ix, []) for anc_ix in ix_ancestor(ix)), set())

        if ap:
            partcount[frozenset(ap)] += 1
            print ix
            print 'full', ' '.join(sorted(af))
            print 'part', ' '.join(sorted(ap))
            print 

            #for aaa in ap:
            #    ageo = admin_partial[ix][aaa]
            #    ageo.intersection(
        else:
            if af:
                fullonlycount += 1
            else:
                nullius.add(ix)

    for k, v in sorted(partcount.items(), key=lambda (k, v): v, reverse=True):
        print ' '.join(sorted(k)), v
    for nn in sorted(nullius):
        c = reduce(lambda a, b: list(a.children())[int(b)], nn, Quad.root()).poly.exterior.coords[0]
        print nn, c[1], c[0]

    print fullonlycount

if __name__ == "__main__":

    import pickle

    #coast_ix = load_coastline(COAST_DATA)
    with open('/home/drew/tmp/coastix') as f:
        coast_ix = pickle.load(f)

    #admin_ix = build_admin_index(ADMIN_DATA)
    with open('/home/drew/tmp/adminix') as f:
        admin_ix = pickle.load(f)

