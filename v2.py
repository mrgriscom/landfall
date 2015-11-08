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
    output = {}
    for ix, geom in coast.iteritems():
        if ix[0] in '23':
            continue

        lines = list(itertools.chain(*map(explode_lines, geom)))
        coverage = dict((ln, set()) for ln in lines)

        # partial admin areas
        for partial_area, admin_bound in admin_partial.get(ix, {}).iteritems():
            assert len(admin_bound) == 1
            admin_bound = admin_bound[0]

            for line in list(coverage.keys()):
                admin_line = admin_bound.intersection(line)
                if admin_line.is_empty or not explode_lines(admin_line):
                    continue

                if line.equals(admin_line):
                    coverage[line].add(partial_area)
                else:
                    remainder = line.difference(admin_bound)
                    cov = coverage.pop(line)
                    for x in explode_lines(remainder):
                        coverage[x] = set(cov)
                    for x in explode_lines(admin_line):
                        coverage[x] = set(cov)
                        coverage[x].add(partial_area)

        # full areas
        af = reduce(lambda a, b: a.union(b), (admin_full.get(anc_ix, []) for anc_ix in ix_ancestor(ix)), set())
        for full_area in af:
            for areas in coverage.values():
                areas.add(full_area)

        output[ix] = coverage
    return output

def verify(output):
    by_admin = collections.defaultdict(list)

    for coverage in output.values():
        for line, areas in coverage.iteritems():
            for aa in areas:
                by_admin[aa].append(line)
            if not areas:
                by_admin['XX'].append(line)

    import csv
    with open('/tmp/lfout.csv', 'w') as f:
        w = csv.DictWriter(f, ['WKT', 'CODE'])
        w.writeheader()
        for area, lines in list(by_admin.iteritems())[:5]:
            w.writerow({'WKT': MultiLineString(lines).wkt, 'CODE': area})

if __name__ == "__main__":

    import pickle

    #coast_ix = load_coastline(COAST_DATA)
    with open('/home/drew/tmp/coastix') as f:
        coast_ix = pickle.load(f)

    #admin_ix = build_admin_index(ADMIN_DATA)
    with open('/home/drew/tmp/adminix') as f:
        admin_ix = pickle.load(f)

    processed = analyze(coast_ix, admin_ix)
    with open('/home/drew/tmp/lf_processed', 'w') as f:
        pickle.dump(processed, f)
