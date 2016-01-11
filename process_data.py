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
import geodesy
import shutil
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--refresh_coast', action='store_true', dest='refresh_coast',
                  help='refresh coast data from openstreetmap')
parser.add_option('--update_coast_corrections', action='store_true', dest='update_coast_corrections',
                  help='reprocess corrections to coastline')
parser.add_option('--update_admin_areas', action='store_true', dest='update_admin_areas',
                  help='reprocess admin areas')

CC_UNCLAIMED = 'X0'
CC_DISPUTED = 'XX'
CC_IGNORE = 'XI'

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

def build_admin_index(path):
    admins = load_geometry(path, 'CODE')
    print '%d areas' % len(admins)

    for k, v in config.expected_unclaimed.iteritems():
        admins['%s-%s' % (CC_UNCLAIMED, k)] = [shapely.wkt.loads(v['bounds'])]
    for k, v in config.incompletely_subdivided.iteritems():
        if v != 'all':
            admins['%s-%s' % (CC_IGNORE, k)] = [shapely.wkt.loads(v)]
    for k, v in config.disputed_areas.iteritems():
        admins['%s-%s' % (CC_DISPUTED, k)] = [shapely.wkt.loads(v['bounds'])]
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

def load_coastline(path, patch_path):
    coastpolys = load_geometry(path)
    patches = load_geometry(patch_path)
    coastpolys = list(apply_patches(coastpolys, patches))
    coastline = list(all_rings(coastpolys))
    full, partial = spatial_decompose(coastline, MAX_DEPTH, verbose=True)
    assert not full
    return partial

def apply_patches(coastpolys, patches):
    for i, poly in enumerate(coastpolys):
        for patch in patches:
            if poly.intersects(patch):
                print 'patching %d' % i
                poly = poly.difference(patch)
        for e in explode_poly(poly):
            yield e

def ix_ancestor(ix):
    for i in xrange(len(ix), -1, -1):
        yield ix[:i]

def is_degenerate_line(line):
    return (line.is_empty or line.length < 1e-9)

def remove_line(ln):
    if is_degenerate_line(ln):
        return True
    polar_cutoff = 90. - geodesy.EPSILON
    _, latmin, _, latmax = ln.bounds
    if latmax < -polar_cutoff or latmin > polar_cutoff:
        # these 'lines' (really just polar singularities) provide no value
        # and slow things down later
        return True
    return False

def analyze(coast, admin):
    admin_full, admin_partial = admin
    output = {}
    for ix, geom in coast.iteritems():
        lines = list(itertools.chain(*map(explode_lines, geom)))
        coverage = dict((ln, set()) for ln in lines if not remove_line(ln))
        if not coverage:
            continue

        # partial admin areas
        for partial_area, admin_bound in admin_partial.get(ix, {}).iteritems():
            for line in list(coverage.keys()):
                admin_line = admin_bound.intersection(line)
                if is_degenerate_line(admin_line):
                    continue

                if line.equals(admin_line):
                    remainder = None
                else:
                    remainder = line.difference(admin_bound)
                    if is_degenerate_line(remainder):
                        remainder = None

                if remainder is None:
                    coverage[line].add(partial_area)
                else:
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

class TopologyError(object):
    pass

class UnclaimedLand(TopologyError):
    def msg(self):
        return 'unexpected unclaimed land'

class DisputedLand(TopologyError):
    def __init__(self, parties):
        self.parties = parties

    def msg(self):
        return 'unexpected disputed area: %s' % ', '.join(sorted(self.parties))

class ParentWithoutSubdivision(TopologyError):
    def __init__(self, parent):
        self.parent = parent

    def msg(self):
        return 'parent w/o subdivision: %s' % self.parent

class SubdivisionWithoutParent(TopologyError):
    def __init__(self, subdiv):
        self.subdiv = subdiv

    def msg(self):
        return 'subdivision not covered by parent: %s' % self.subdiv

def verify(output, admin_info, verbose=False):
    subdiv_parents = set(e['parent'] for e in admin_info.values() if e['type'] == 'subdivision')

    for ix, coverage in sorted(output.iteritems()):
        issues = set()
        problem_lines = set()
        for line, areas in coverage.iteritems():
            special = {}
            for spec in (CC_UNCLAIMED, CC_IGNORE, CC_DISPUTED):
                special[spec] = set(a for a in areas if a.startswith('%s-' % spec))
                areas -= special[spec]

            # check for issues

            if not areas and not special[CC_UNCLAIMED]:
                issues.add(UnclaimedLand())

            parents = set(admin_info[a]['parent'] for a in areas if admin_info[a]['type'] == 'subdivision')
            for a in areas:
                info = admin_info[a]

                if info['type'] == 'subdivision' and info['parent'] not in areas:
                    issues.add(SubdivisionWithoutParent(a))

                if a in subdiv_parents and a not in parents:
                    cci = '%s-%s' % (CC_IGNORE, a)
                    if cci in special[CC_IGNORE] or config.incompletely_subdivided.get(a) == 'all':
                        pass
                    else:
                        issues.add(ParentWithoutSubdivision(a))
                        problem_lines.add(line)
            areas -= parents

            if len(areas) > 1:
                valid_dispute = None

                if len(areas) == 2:
                    for defacto, recognized in config.defacto.iteritems():
                        if areas.issubset(set([defacto] + recognized)):
                            valid_dispute = defacto
                            break

                if not valid_dispute:
                    for d in special[CC_DISPUTED]:
                        did = d.split('-')[-1]
                        if set(config.disputed_areas[did]['parties']) == areas:
                            valid_dispute = d
                            break

                if valid_dispute:
                    # need to modify original 'areas' in the coverage dict, so can't set to new set
                    areas -= areas
                    areas.add(valid_dispute)
                else:
                    issues.add(DisputedLand(areas))
                    problem_lines.add(line)

        # print warnings
        if issues:
            coords = tuple(reversed(Quad.from_ix(ix).poly.centroid.coords[0]))

            if verbose and problem_lines:
                print '-----'
                for ln in problem_lines:
                    print ln

            for iss in issues:
                print '%s; near %.3f %.3f' % (iss.msg(), coords[0], coords[1])

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

            # quick fix for antarctica -- isn't converted to wkt correctly
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

def process(options):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    tmp_dir = os.path.join(data_dir, 'tmp')

    admin_info_path = os.path.join(data_dir, 'admin_index.csv')
    admin_dir = os.path.join(data_dir, 'admin_areas')
    noland_dir = os.path.join(data_dir, 'no_land')

    coast_dir = os.path.join(tmp_dir, 'land-polygons-complete-4326')
    coast_csv = os.path.join(tmp_dir, 'coast.csv')
    noland_csv = os.path.join(tmp_dir, 'noland.csv')
    admin_csv = os.path.join(tmp_dir, 'admin.csv')
    admin_ix_path = os.path.join(tmp_dir, 'admin_ix')
    coast_ix_path = os.path.join(tmp_dir, 'coast_ix')
    final_output_path = os.path.join(tmp_dir, 'tagged_coastline')

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    wipeouts = set()
    if options.refresh_coast:
        wipeouts.update([coast_dir, coast_csv, coast_ix_path, final_output_path])
    if options.update_coast_corrections:
        wipeouts.update([noland_csv, coast_ix_path, final_output_path])
    if options.update_admin_areas:
        wipeouts.update([admin_csv, admin_ix_path, final_output_path])
    for w in wipeouts:
        if os.path.exists(w):
            print 'wiping out %s' % w
            if w in (coast_dir,):
                shutil.rmtree(w)
            else:
                os.remove(w)

    if not os.path.exists(coast_dir):
        print 'Downloading coastline data...'
        zip_path = tempfile.mktemp()
        os.popen('wget -O %s http://data.openstreetmapdata.com/land-polygons-complete-4326.zip' % zip_path)
        os.popen('unzip -d %s %s' % (tmp_dir, zip_path))
    else:
        print 'Reusing existing coastline data'

    if not os.path.exists(coast_csv):
        print 'Extracting coastline data...'
        shp_to_csv(find_with_ext(coast_dir, 'shp'), coast_csv)
    else:
        print 'Reusing extracted coastline data'

    if not os.path.exists(noland_csv):
        print 'Extracting coastline patch data...'
        shp_to_csv(find_with_ext(noland_dir, 'shp'), noland_csv)
    else:
        print 'Reusing extracted coastline patch data'

    if not os.path.exists(admin_csv):
        print 'Extracting admin data...'
        shp_to_csv(find_with_ext(admin_dir, 'shp'), admin_csv)
    else:
        print 'Reusing extracted admin data'

    if not os.path.exists(admin_ix_path):
        print 'Building admin area index...'
        admin_ix = build_admin_index(admin_csv)
        with open(admin_ix_path, 'w') as f:
            pickle.dump(admin_ix, f)
    else:
        print 'Loading admin area index...'
        with open(admin_ix_path) as f:
            admin_ix = pickle.load(f)

    if not os.path.exists(coast_ix_path):
        print 'Building coastline index...'
        coast_ix = load_coastline(coast_csv, noland_csv)
        with open(coast_ix_path, 'w') as f:
            pickle.dump(coast_ix, f)
    else:
        print 'Loading coastline index...'
        with open(coast_ix_path) as f:
            coast_ix = pickle.load(f)

    print 'Tagging coastline with admin regions...'
    processed = analyze(coast_ix, admin_ix)
    print 'Checking result for problems...'
    admin_info = load_admin_info(admin_info_path)
    verify(processed, admin_info)
    print 'Check complete; writing final output...'
    with open(final_output_path, 'w') as f:
        pickle.dump(processed, f)

def load_admin_info(path):
    info = {}
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            info[row['code']] = row
    return info

if __name__ == "__main__":

    (options, args) = parser.parse_args()
    process(options)



