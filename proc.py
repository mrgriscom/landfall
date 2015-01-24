import shapely
import shapely.wkt
import shapely.ops
import shapely.validation
from shapely.geometry import Polygon, MultiPolygon
import sys
import csv
import itertools
import collections

Extent = collections.namedtuple('Extent', 'x0 x1 y0 y1')

def display(poly):
    from matplotlib import pyplot
    from matplotlib.patches import Circle
    from shapely.geometry import Polygon
    from descartes.patch import PolygonPatch

    def plot_coords(ax, ob):
        x, y = ob.xy
        ax.plot(x, y, 'o', color='#999999', zorder=1)
    
    fig = pyplot.figure(1, figsize=[5, 5], dpi=90)

    # 3: invalid polygon, ring touch along a line
    ax = fig.add_subplot(111)

    polys = poly if isinstance(poly, MultiPolygon) else [poly]

    for poly in polys:
        for inter in poly.interiors:
            plot_coords(ax, inter)
        plot_coords(ax, poly.exterior)

        patch = PolygonPatch(poly, facecolor='#888888', edgecolor='#222222')
        ax.add_patch(patch)

    #ax.set_aspect(1)
    pyplot.show()
    

def load():
    csv.field_size_limit(sys.maxsize)
    data = []
    with open('/home/drew/tmp/land_polygons.csv') as f:
        r = csv.DictReader(f)
        for row in r:
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
            data.append(poly)
    return data

def mkquad(ext):
    return Polygon([(ext.x0, ext.y0), (ext.x1, ext.y0), (ext.x1, ext.y1), (ext.x0, ext.y1)])

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

def test(data, quad, merge=False):
    if merge and len(data) > 1:
        data = [shapely.ops.cascaded_union(data)]

    subdata = []
    for poly in data:
        if quad.intersects(poly):
            if quad.within(poly):
                return True
            else:
                subdata.append(quad.intersection(poly))
    return subdata if subdata else False


def quadtree(ext, maxdepth, data, consolidate_final=True, depth=0):
    print depth, ext.x0, ext.y0, ext.x1-ext.x0, ext.y1-ext.y0, len(data)

    result = test(data, mkquad(ext), consolidate_final and depth == maxdepth)
    if result is False:
        return -1
    elif result is True:
        return 1
    else:
        if depth == maxdepth:
            return 0
        else:
            children = [quadtree(child, maxdepth, result, consolidate_final, depth + 1) for child in quadchildren(ext)]
            if children == [1, 1, 1, 1]:
                # can happen because polygons are only consolidated once we hit max_depth
                children = 1
            return children


def proc_index(node, handler, depth=0, tile=(0, 0)):
    if hasattr(node, '__iter__'):
        for child, subtile in zip(node, [(2*tile[0] + i, 2*tile[1] + j) for j in xrange(2) for i in xrange(2)]):
            proc_index(child, handler, depth + 1, subtile)
    else:
        handler(node, depth, tile)

def run(data):
    DEPTH = 15
    EXTENT = Extent(-180, 180, 90, -270)
    #EXTENT = Extent(0, 45, 10, -35)

    return quadtree(EXTENT, DEPTH, data, DEPTH > 8)

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


if __name__ == "__main__":
    data = load()
    print 'loaded'
    run(data)
