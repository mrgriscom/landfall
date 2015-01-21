import shapely
import shapely.wkt
import shapely.ops
from shapely.geometry import Polygon
import sys
import csv
import itertools

def load():
    csv.field_size_limit(sys.maxsize)
    with open('/tmp/coast2/blah/land_polygons.csv') as f:
        r = csv.DictReader(f)
        return [shapely.wkt.loads(row['WKT']) for row in r]

def mkquad(x0, x1, y0, y1):
    return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])

def pr(s):
    sys.stdout.write(s)
    sys.stdout.flush()

def test(data, quad, merge=False):
    if merge:
        data = [shapely.ops.cascaded_union(data)]

    subdata = []
    for poly in data:
        if quad.intersects(poly):
            if quad.within(poly):
                return True
            else:
                subdata.append(quad.intersection(poly))
    return subdata if subdata else False


def quadtree(x0, x1, y0, y1, maxdepth, data, depth=0, tile=(0, 0)):
    print '... %s %s %s' % (depth, tile, len(data))
    print x0, x1, y0, y1
    result = test(data, mkquad(x0, x1, y0, y1), depth == maxdepth)
    if result is False:
        return
    elif result is True:
        yield ((depth, tile[0], tile[1]), True)
    else:
        if depth == maxdepth:
            yield ((depth, tile[0], tile[1]), False)
        else:
            xmid = .5*(x0+x1)
            ymid = .5*(y0+y1)
            for e in itertools.chain(
                quadtree(x0, xmid, y0, ymid, maxdepth, result, depth + 1, (2*tile[0]+0, 2*tile[1]+0)),
                quadtree(xmid, x1, y0, ymid, maxdepth, result, depth + 1, (2*tile[0]+1, 2*tile[1]+0)),
                quadtree(x0, xmid, ymid, y1, maxdepth, result, depth + 1, (2*tile[0]+0, 2*tile[1]+1)),
                quadtree(xmid, x1, ymid, y1, maxdepth, result, depth + 1, (2*tile[0]+1, 2*tile[1]+1)),
            ):
                yield e


def run(data):
    DEPTH = 7
    DIM = 2**DEPTH
    BITMAP = [[0 for i in xrange(DIM)] for j in xrange(DIM)]
    #EXTENT = [-180, 0, -90, 90]
    EXTENT = [-80, -55, 30, 55]

    for e in quadtree(EXTENT[0], EXTENT[1], EXTENT[2], EXTENT[3], DEPTH, data):
        (z, x, y), st = e
        celldim = 2**(DEPTH - z)
        for i in xrange(celldim):
            for j in xrange(celldim):
                BITMAP[DIM-1-(y*celldim + j)][x*celldim + i] = 2 if st else 1

    print '\n'.join(''.join('.xO'[col] for col in row) for row in BITMAP)


if __name__ == "__main__":
    data = load()
    print 'loaded'
    run(data)
