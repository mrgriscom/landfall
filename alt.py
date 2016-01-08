from deprecated_landfall import *
import math
import geodesy
from shapely.geometry import LineString, LinearRing

def all_rings(polys):
    for p in polys:
        for q in explode_poly(p):
            yield q.exterior
            for r in q.interiors:
                yield r

def all_vertices(polys):
    for r in all_rings(polys):
        for c in r.coords:
            yield c

def test(p0, polys):
    ctx = distbear_init(p0)
    for i, v in enumerate(all_vertices(polys)):
        dist, bear = distbear(ctx, swap(v))
        print dist, bear
        if i % 10000 == 0:
            print i

def distbear_init(p0):
    v0 = geodesy.ll_to_ecefu(p0)
    vnorth, veast = geodesy.orientate(v0)
    return v0, vnorth, veast

def distbear(ctx, p1):
    v0, vnorth, veast = ctx
    v1 = geodesy.ll_to_ecefu(p1)

    (vo, kcos) = geodesy.vortho(v0, v1)
    ksin = geodesy.vlen(vo)
    dist = geodesy.EARTH_MEAN_RAD * math.atan2(ksin, kcos)

    if ksin < 1e-9:
        #antipodal
        bear = None
    else:
        bear = geodesy._xy_to_bearing(geodesy.dotp(vo, veast), geodesy.dotp(vo, vnorth))

    return dist, bear

EARTH_CIRCUMF = 2*math.pi*geodesy.EARTH_MEAN_RAD

# TODO: bearing ranges crossing north
def postings(p0, data, bearing_res, bearing0, bearing1, near_dist):
    pixels = int(round((bearing1 - bearing0) / bearing_res))
    output = [(-1, set())] * pixels

    def px(bear):
        return int((bear - bearing0) / bearing_res)

    def update(dist, px, areas):
        if dist < near_dist:
            return
        if px < 0 or px >= pixels:
            return
        if output[px][0] < 0 or dist < output[px][0]:
            output[px] = (dist, areas)

    def fill((dist, bear, x), (prev_dist, prev_bear, prev_x), areas):
        if abs(x - prev_x) <= 1 or abs(x - prev_x) >= pixels - 1:
            return

        minx = min(x, prev_x)
        maxx = max(x, prev_x)
        if maxx - minx > pixels / 2:
            minx, maxx = maxx, minx + pixels

        for i in xrange(minx + 1, maxx):
            # add some law of cosines shit to vary dist
            update(dist, i % pixels, areas)

    def segments():
        for seg in data.values():
            for line, areas in seg.iteritems():
                yield (line, areas)

    i = 0
    ctx = distbear_init(p0)
    for seg, areas in segments():
        is_ring = isinstance(seg, LinearRing)
        coords = seg.coords
        #if is_ring:
        #    coords = list(seg.coords)[:-1]
        #else:
        #    coords = seg.coords

        prev = None
        anti_prev = None
        first = None
        anti_first = None
        for v in coords:
            dist, bear = distbear(ctx, swap(v))
            bear = bear % 360.
            x = px(bear)
            anti_dist = EARTH_CIRCUMF - dist
            anti_bear = (bear + 180.) % 360.
            anti_x = px(anti_bear)

            update(dist, x, areas)
            update(anti_dist, anti_x, areas)

            if prev is not None:
                fill((dist, bear, x), prev, areas)
                fill((anti_dist, anti_bear, anti_x), anti_prev, areas)
            else:
                first = dist, bear, x
                anti_first = anti_dist, anti_bear, anti_x

            prev = dist, bear, x
            anti_prev = anti_dist, anti_bear, anti_x

            i += 1
            if i % 10000 == 0:
                print i

        #if is_ring:
        #    fill(first, prev)
        #    fill(anti_first, anti_prev)

    return output



def enclosing_circle_for_bbox(ll, ur):
    """Choose a circle that encloses the lat/lon bounding box. Note that this will
    not always choose the smallest enclosing circle, and there are certain
    pathological cases that will yield a circle much larger than the bounding box
    (such as high polar bboxes with longitude spans approaching 180 degrees), but
    in practice with boxes of roughly equal latitude/longitude span it's good
    enough."""
    latmin, lonmin = ll
    latmax, lonmax = ur
    assert latmin <= latmax
    assert (lonmax - lonmin) % 360. <= 180.

    vA = geodesy.ll_to_ecefu((latmin, lonmin))
    vB = geodesy.ll_to_ecefu((latmin, lonmax))
    vC = geodesy.ll_to_ecefu((latmax, lonmin))
    vD = geodesy.ll_to_ecefu((latmax, lonmax))
    vAB = geodesy.vdiff(vB, vA)
    vAC = geodesy.vdiff(vC, vA)
    vCD = geodesy.vdiff(vD, vC)
    
    degenerate_points = None
    if geodesy.vlen(vAC) < geodesy.EPSILON:
        # d-lat = 0
        degenerate_points = [vA, vB]
    elif geodesy.vlen(vAB) < geodesy.EPSILON:
        # south pole or d-lon = 0
        if geodesy.vlen(vCD) < geodesy.EPSILON:
            # north pole or d-lon = 0
            degenerate_points = [vA, vC]
        else:
            vects = [vCD, vAC]
    else:
        vects = [vAB, vAC]

    if degenerate_points:
        dp0, dp1 = degenerate_points
        mean = geodesy.vscale(geodesy.vadd(dp0, dp1), .5)
        if geodesy.vlen(mean) < geodesy.EPSILON:
            # degenerate points are antipodal
            ortho = geodesy.ll_to_ecefu((.5*(latmin+latmax), .5*(lonmin+lonmax)))
        else:
            ortho = geodesy.vnorm(mean)
    else:
        ortho = geodesy.vnorm(geodesy.crossp(*vects))
    center = geodesy.ecefu_to_ll(ortho)
    dist = geodesy.distance(center, (latmin, lonmin))
    return center, dist

def bounds_to_circle(p0, center, radius):
    dist = geodesy.distance(p0, center)
    mindist = max(0, dist - radius)
    maxdist = min(.5 * EARTH_CIRCUMF, dist + radius)

    sin_ratio = math.sin(radius / geodesy.EARTH_MEAN_RAD) / math.sin(dist / geodesy.EARTH_MEAN_RAD)
    if abs(sin_ratio) < 1:
        angle = math.degrees(math.asin(sin_ratio))
        bearing = geodesy.bearing(p0, center)
        minbear = geodesy.anglenorm(bearing - angle)
        maxbear = geodesy.anglenorm(bearing + angle)
    else:
        minbear = -180.
        maxbear = 180.

    return mindist, maxdist, minbear, maxbear
