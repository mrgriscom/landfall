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




def lat_for_isodist_tangent(p0, lon):
    dlon = geodesy.anglenorm(lon - p0[1])
    x, y, z = geodesy.ll_to_ecefu((p0[0], dlon))
    if abs(x) < geodesy.EPSILON:
        return None
    if x < 0:
        x = -x
        z = -z
    return math.degrees(math.atan2(z, x))

def lon_for_greatcircle_maxlat(p0, lat):
    z = math.sin(math.radians(p0[0]))
    if abs(lat) > 90. - geodesy.EPSILON:
        x = 0
    else:
        x = z / math.tan(math.radians(lat))
    y2 = 1. - x**2 - z**2
    if y2 < 0:
        return []
    y = y2**.5
    dlon = math.degrees(math.atan2(y, x))
    return [geodesy.anglenorm(p0[1] + k*dlon) for k in (-1, 1)]

# inclusive
def in_lon_range(lon, lonmin, lonmax):
    dlon = (lon - lonmin) % 360.
    return dlon >= 0 and dlon <= (lonmax - lonmin) % 360.

DEBUG = True
def log(s):
    if DEBUG:
        print s

def minmax_bbox_dist(p0, ll, ur):
    lat, lon = p0
    latmin, lonmin = ll
    latmax, lonmax = ur

    anti_lon = lon + 180.
    anti_latmin = -latmax
    anti_latmax = -latmin
    anti_lonmin = lonmin + 180.
    anti_lonmax = lonmax + 180.

    def key_points():
        for latbound in (latmin, latmax):
            for lonbound in (lonmin, lonmax):
                yield (latbound, lonbound)
        for _lon in (lon, anti_lon):
            if in_lon_range(_lon, lonmin, lonmax):
                for latbound in (latmin, latmax):
                    yield (latbound, _lon)
        for _lon in (lonmin, lonmax):
            key_lat = lat_for_isodist_tangent(p0, _lon)
            if key_lat is not None and key_lat > latmin and key_lat < latmax:
                yield (key_lat, _lon)

    log(list(key_points()))
    dists = [geodesy.distance(p0, key_point) for key_point in key_points()]

    if lat >= latmin and lat <= latmax and in_lon_range(lon, lonmin, lonmax):
        log('min at pode')
        mindist = 0
    else:
        log('min at %s' % str(min(key_points(), key=lambda e: geodesy.distance(p0, e))))
        mindist = min(dists)

    if lat >= anti_latmin and lat <= anti_latmax and in_lon_range(lon, anti_lonmin, anti_lonmax):
        log('max at antipode')
        maxdist = .5 * EARTH_CIRCUMF
    else:
        log('max at %s' % str(max(key_points(), key=lambda e: geodesy.distance(p0, e))))
        maxdist = max(dists)

    return mindist, maxdist

def minmax_bbox_bearing(p0, ll, ur):
    lat, lon = p0
    latmin, lonmin = ll
    latmax, lonmax = ur

    anti_latmin = -latmax
    anti_latmax = -latmin
    anti_lonmin = lonmin + 180.
    anti_lonmax = lonmax + 180.
    
    # handles all boundary cases
    if ((lat >= latmin and lat <= latmax and in_lon_range(lon, lonmin, lonmax)) or
        (lat >= anti_latmin and lat <= anti_latmax and in_lon_range(lon, anti_lonmin, anti_lonmax))):
        log('(anti)pode')
        return None

    def key_points():
        for latbound, vtype in zip((latmin, latmax), ('b', 't')):
            for lonbound, htype in zip((lonmin, lonmax), ('l', 'r')):
                yield ((latbound, lonbound), vtype + htype)
        key_lat = None
        if latmin > 0:
            key_lat = latmin
        elif latmax < 0:
            key_lat = latmax
        if key_lat is not None:
            for key_lon, type in zip(lon_for_greatcircle_maxlat(p0, key_lat), ('tang_l', 'tang_r')):
                if in_lon_range(key_lon, lonmin, lonmax):
                    yield ((key_lat, key_lon), type)

    def inflection_points():
        for key_point, type in key_points():
            if type in ('tang_l', 'tang_r'):
                yield key_point
            else:
                rev_bearing = geodesy.bearing(key_point, p0)
                quadrant = math.floor(rev_bearing / 90.) % 4
                if quadrant % 2 == {'tl': 0, 'tr': 1, 'bl': 1, 'br': 0}[type]:
                    yield key_point

    for infl_point in inflection_points():
        bearing = geodesy.bearing(p0, infl_point)
        print infl_point, bearing

    return


    if lat > latmin and lat < latmax and in_lon_range(lon, lonmin, lonmax):
        log('min at pode')
        mindist = 0
    else:
        log('min at %s' % str(min(key_points(), key=lambda e: geodesy.distance(p0, e))))
        mindist = min(dists)

    if lat > anti_latmin and lat < anti_latmax and in_lon_range(lon, anti_lonmin, anti_lonmax):
        log('max at antipode')
        maxdist = .5 * EARTH_CIRCUMF
    else:
        log('max at %s' % str(max(key_points(), key=lambda e: geodesy.distance(p0, e))))
        maxdist = max(dists)

    return mindist, maxdist

