import shapely
import shapely.wkt
import shapely.ops
import shapely.validation
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection, box
import sys
import csv
import itertools
import collections
import tempfile
import json

Extent = collections.namedtuple('Extent', 'x0 x1 y0 y1')

def explode_poly(poly):
    if isinstance(poly, Polygon):
        return [poly]
    elif isinstance(poly, MultiPolygon):
        return poly
    elif isinstance(poly, GeometryCollection):
        return [e for e in poly if isinstance(e, Polygon)]
    else:
        return []

def display(poly):
    from matplotlib import pyplot
    from matplotlib.patches import Circle
    from shapely.geometry import Polygon
    from descartes.patch import PolygonPatch

    def plot_coords(ax, ob):
        x, y = ob.xy
        ax.plot(x, y, 'o', color='#999999', zorder=1)
    
    fig = pyplot.figure(1, figsize=[5, 5], dpi=90)

    ax = fig.add_subplot(111)

    for poly in explode_poly(poly):
        for inter in poly.interiors:
            plot_coords(ax, inter)
        plot_coords(ax, poly.exterior)

        patch = PolygonPatch(poly, facecolor='#888888', edgecolor='#222222')
        ax.add_patch(patch)

    #ax.set_aspect(1)
    pyplot.show()

COAST_DATA = '/tmp/coastcsv/land_polygons.csv'
#COAST_DATA = '/home/drew/tmp/land_polygons.csv'
#ADMIN_DATA = '/home/drew/tmp/admin_polygons.csv'
def load(path, separate=False):
    csv.field_size_limit(sys.maxsize)
    data = {} if separate else []
    with open(path) as f:
        r = csv.DictReader(f)
        for i, row in enumerate(r):
            poly = shapely.wkt.loads(row['WKT'])
            if not poly.is_valid:
                #display(poly)
                reason = shapely.validation.explain_validity(poly).lower()
                print reason
                if 'self-intersection' in reason:
                    poly = poly.buffer(0.)
                elif 'too few points' in reason:
                    pass
                else:
                    assert False, poly
            set_complexity(poly)
            if separate:
                data[i] = [poly]
            else:
                data.append(poly)
    return data

def mkquad(ext):
    return box(ext.x0, ext.y0, ext.x1, ext.y1)

ROOT = Extent(-180, 180, 90, -270)

def quadchildren(ext):
    xmid = .5 * (ext.x0 + ext.x1)
    ymid = .5 * (ext.y0 + ext.y1)
    yield Extent(ext.x0, xmid, ext.y0, ymid)
    yield Extent(xmid, ext.x1, ext.y0, ymid)
    yield Extent(ext.x0, xmid, ymid, ext.y1)
    yield Extent(xmid, ext.x1, ymid, ext.y1)

def pr(s):
    sys.stdout.write(s)
    sys.stdout.flush()

def poly_complexity(poly):
    cmplx = lambda ring: len(ring.coords)
    num_vertices = lambda p: cmplx(p.exterior) + sum(map(cmplx, p.interiors))
    if isinstance(poly, MultiPolygon):
        return sum(map(num_vertices, poly))
    else:
        return num_vertices(poly)

def set_complexity(poly):
    poly.complexity = sum(map(poly_complexity, explode_poly(poly)))

SUBDIV_CUTOFF = 10
def test(data, quad):
    subdata = []
    for poly in data:
        if quad.intersects(poly):
            if quad.within(poly):
                return True
            else:
                if poly.complexity > SUBDIV_CUTOFF and not quad.contains(poly):
                    subpoly = quad.intersection(poly)
                    set_complexity(subpoly)
                else:
                    subpoly = poly
                subdata.append(subpoly)
    return subdata if subdata else False

COVERAGE_FULL = 1
COVERAGE_PARTIAL = 0
COVERAGE_NONE = -1

class QTNode(object):
    #__slots__ = ['children', 'land', 'aux_partial', 'aux_full']
    def __init__(self):
        self.children = None
        self.land = None
        self.aux_partial = []
        self.aux_full = []

    def tidy_up(self):
        for attr in ('aux_full', 'aux_partial'):
            val = getattr(self, attr) or []
            val = tuple(val) if val else None
            setattr(self, attr, val)

def quadtree(ext, maxdepth, data, depth=0, ix='', cache_data=lambda ix, ext, data: None):
    print (ix or '~'), len(data)
    ret = QTNode()

    quad = mkquad(ext)
    result = test(data, quad)
    if result is False:
        # current extent is all water
        return None

    if result is True:
        ret.land = COVERAGE_FULL
    else:
        if depth == maxdepth:
            ret.land = COVERAGE_PARTIAL
            cache_data(ix, ext, result)
        else:
            quadtree_descend(ret, maxdepth, ext, result, ix, cache_data)

    ret.tidy_up()
    return ret

def quadtree_descend(node, maxdepth, ext, data, ix, cache_data):
    node.land = None
    depth = len(ix)
    node.children = [quadtree(child, maxdepth, data, depth + 1, ix + s, cache_data) for child, s in zip(quadchildren(ext), '0123')]

def dump_ix(root, path=None):
    import pickle
    if not path:
        path = tempfile.mktemp()
    with open(path, 'w') as f:
        pickle.dump(root, f)
    return path

def load_ix(path):
    import pickle
    with open(path) as f:
        return pickle.load(f)

def proc_index(node, handler, depth=0, tile=(0, 0)):
    if not node:
        return
    if node.children:
        for child, subtile in zip(node.children, [(2*tile[0] + i, 2*tile[1] + j) for j in xrange(2) for i in xrange(2)]):
            proc_index(child, handler, depth + 1, subtile)
    else:
        handler(node, depth, tile)

def build_index(raw, extent, max_depth, cache_data=None):
    return quadtree(extent, max_depth, raw, cache_data=cache_data)

GLOBAL = Extent(-180, 180, 90, -270)

def render(ix, depth):
    dim = 2**depth
    BITMAP = [[0 for i in xrange(dim)] for j in xrange(dim)]
    def paint(val, z, (x, y)):
        c = {None: 0, 0: 128, 1: 255}[val.land]
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

if __name__ == "__main__":
    max_depth = int(sys.argv[1])

    coast = load(COAST_DATA)
    #admin = load(ADMIN_DATA, True)
    print 'loaded'
    ix = build_index(coast, GLOBAL, max_depth)
    print 'index built'
    print dump_ix(ix)
