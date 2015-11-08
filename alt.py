from landfall import *
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
    output = [None] * pixels

    def px(bear):
        return int((bear - bearing0) / bearing_res)

    def update(dist, px, areas):
        if dist < near_dist:
            return
        if px < 0 or px >= pixels:
            return
        if output[px] is None or dist < output[px][0]:
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
