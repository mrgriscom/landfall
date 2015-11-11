import pickle
import os
import os.path
import tempfile
import csv
import sys
import collections
import itertools
import shapely
import shapely.wkt
import shapely.ops
import shapely.validation
import shapely.geos
from shapely.geometry import LineString, LinearRing, MultiLineString, Polygon, MultiPolygon, GeometryCollection, box
import config

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

    @staticmethod
    def from_ix(ix):
        return reduce(lambda a, b: list(a.children())[int(b)], ix, Quad.root())

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
        if ix[0:1] in ('2', '3'):
            return

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

def explode_poly(poly):
    if isinstance(poly, Polygon):
        return [poly]
    elif isinstance(poly, MultiPolygon):
        return poly
    elif isinstance(poly, GeometryCollection):
        return [e for e in poly if isinstance(e, Polygon)]
    else:
        return []

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

    for k, v in config.expected_unclaimed.iteritems():
        admins['%s-%s' % (config.cc_unclaimed, k)] = [shapely.wkt.loads(v['bounds'])]
    for k, v in config.incompletely_subdivided.iteritems():
        if v != 'all':
            admins['%s-%s' % (config.cc_ignore, k)] = [shapely.wkt.loads(v)]
    for k, v in config.disputed_areas.iteritems():
        admins['%s-%s' % (config.cc_disputed, k)] = [shapely.wkt.loads(v['bounds'])]
    print '%d areas (with supplemental)' % len(admins)

    ixfull = collections.defaultdict(set)
    ixpartial = collections.defaultdict(dict)

    for k, v in sorted(admins.items()):
        full, partial = spatial_decompose(v, MAX_DEPTH)
        print k, len(full), len(partial)
        for e in full:
            ixfull[e].add(k)
        for e, f in partial.iteritems():
            assert len(f) == 1
            ixpartial[e][k] = f[0]

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
        lines = list(itertools.chain(*map(explode_lines, geom)))
        coverage = dict((ln, set()) for ln in lines)

        # partial admin areas
        for partial_area, admin_bound in admin_partial.get(ix, {}).iteritems():
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

def verify(output, admin_info):
    subdiv_parents = set(e['parent'] for e in admin_info.values() if e['type'] == 'subdivision')

    for ix, coverage in sorted(output.iteritems()):
        issues = set()
        for line, areas in coverage.iteritems():
            # TODO fix this in indexing phase
            if len(line.coords) == 2 and line.length < 1e-9:
                continue

            special = {}
            for spec in (config.cc_unclaimed, config.cc_ignore, config.cc_disputed):
                special[spec] = set(a for a in areas if a.startswith('%s-' % spec))
                areas -= special[spec]

            # check for issues
            if not areas and not special[config.cc_unclaimed]:
                issues.add((config.cc_unclaimed,))

            parents = set(filter(None, (admin_info[a]['parent'] for a in areas)))
            for a in areas:
                info = admin_info[a]

                if info['type'] == 'subdivision' and info['parent'] not in areas:
                    issues.add(('XPAR', a))

                if a in subdiv_parents and a not in parents:
                    cci = '%s-%s' % (config.cc_ignore, a)
                    if cci in special[config.cc_ignore] or config.incompletely_subdivided.get(a) == 'all':
                        pass
                    else:
                        issues.add(('PPAR', a))
            areas -= parents

            if len(areas) > 1:
                valid_dispute = None
                for d in special[config.cc_disputed]:
                    did = d.split('-')[-1]
                    if set(config.disputed_areas[did]['parties']) == areas:
                        valid_dispute = d
                        break
                if valid_dispute:
                    areas.add(valid_dispute)
                else:
                    issues.add((config.cc_disputed, tuple(sorted(areas))))

        # print warnings
        if issues:
            coords = tuple(reversed(Quad.from_ix(ix).poly.centroid.coords[0]))
            def err(s):
                print '%s near %.3f %.3f' % (s, coords[0], coords[1])

            for i in issues:
                type = i[0]
                if type == config.cc_unclaimed:
                    err('unexpected unclaimed land')
                elif type == 'XPAR':
                    err('subdivision not covered by parent [%s]' % i[1])
                elif type == 'PPAR':
                    err('parent w/o a subdivision [%s]' % i[1])
                elif type == config.cc_disputed:
                    err('disputed area [%s]' % ', '.join(i[1]))


def load_geometry(path, key=None):
    """if coalesce=true, all polys are part of one big multipolygon; if false,
    each represents a separate entity"""
    csv.field_size_limit(sys.maxsize)
    data = [] if key is None else {}
    with open(path) as f:
        r = csv.DictReader(f)
        for i, row in enumerate(r):
            if (i+1) % 10000 == 0:
                print 'loaded %d' % (i+1)

            # quick fix for antarctica
            if row['WKT'] == 'POLYGON ((-180 -90,-180 -60,180 -60,180 -90))':
                row['WKT'] = 'POLYGON ((-180 -90,-180 -60,180 -60,180 -90,-180 -90))'

            try:
                poly = shapely.wkt.loads(row['WKT'])
            except shapely.geos.ReadingError:
                assert False, 'invalid geometry'
            if not poly.is_valid:
                reason = shapely.validation.explain_validity(poly).lower()
                print reason
                if 'self-intersection' in reason:
                    poly = poly.buffer(0.)
                elif 'too few points' in reason:
                    pass
                else:
                    assert False, poly
            if key:
                if row[key] in data:
                    print 'warning: dup key %s' % row[key]
                data[row[key]] = [poly]
            else:
                data.append(poly)
    return data

def find_with_ext(path, ext):
    matches = [f for f in os.listdir(path) if f.endswith('.%s' % ext)]
    assert len(matches) == 1
    return os.path.join(path, matches[0])

def shp_to_csv(shppath, csvpath):
    tmp = tempfile.mktemp()  # mkdtemp actually creates the dir, which we don't want
    os.popen('ogr2ogr -f csv -lco GEOMETRY=AS_WKT -progress %s %s' % (tmp, shppath))
    os.popen('mv "%s" "%s"' % (find_with_ext(tmp, 'csv'), csvpath))

def process():
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    tmp_dir = os.path.join(data_dir, 'tmp')

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    coast_dir = os.path.join(tmp_dir, 'land-polygons-complete-4326')
    if not os.path.exists(coast_dir):
        print 'Downloading coastline data...'
        zip_path = tempfile.mktemp()
        os.popen('wget -O %s http://data.openstreetmapdata.com/land-polygons-complete-4326.zip' % zip_path)
        os.popen('unzip -d %s %s' % (tmp_dir, zip_path))
    else:
        print 'Reusing existing coastline data'

    coast_csv = os.path.join(tmp_dir, 'coast.csv')
    if not os.path.exists(coast_csv):
        print 'Extracting coastline data...'
        shp_to_csv(find_with_ext(coast_dir, 'shp'), coast_csv)
    else:
        print 'Resuing extracted coastline data'

    admin_dir = os.path.join(data_dir, 'admin_areas')
    admin_csv = os.path.join(tmp_dir, 'admin.csv')
    if not os.path.exists(admin_csv):
        print 'Extracting admin data...'
        shp_to_csv(find_with_ext(admin_dir, 'shp'), admin_csv)
    else:
        print 'Resuing extracted admin data'

    admin_ix_path = os.path.join(tmp_dir, 'admin_ix')
    if not os.path.exists(admin_ix_path):
        print 'Building admin area index...'
        admin_ix = build_admin_index(admin_csv)
        with open(admin_ix_path, 'w') as f:
            pickle.dump(admin_ix, f)
    else:
        print 'Loading admin area index...'
        with open(admin_ix_path) as f:
            admin_ix = pickle.load(f)

    coast_ix_path = os.path.join(tmp_dir, 'coast_ix')
    if not os.path.exists(coast_ix_path):
        print 'Building coastline index...'
        coast_ix = load_coastline(coast_csv)
        with open(coast_ix_path, 'w') as f:
            pickle.dump(coast_ix, f)
    else:
        print 'Loading coastline index...'
        with open(coast_ix_path) as f:
            coast_ix = pickle.load(f)

    #final_output_path = os.path.join(tmp_dir, 'tagged_coastline')
    processed = analyze(coast_ix, admin_ix)
    admin_info_path = os.path.join(data_dir, 'admin_index.csv')
    admin_info = load_admin_info(admin_info_path)
    verify(processed, admin_info)
    #with open(final_output_path, 'w') as f:
    #    pickle.dump(processed, f)

    import pdb;pdb.set_trace()

def load_admin_info(path):
    info = {}
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            info[row['code']] = row
    return info

if __name__ == "__main__":

    process()



