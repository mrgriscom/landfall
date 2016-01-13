import math
import geodesy
from shapely.geometry import LineString, LinearRing
import process_data as pd
import Queue
import itertools

EARTH_CIRCUMF = 2*math.pi*geodesy.EARTH_MEAN_RAD
EARTH_FARTHEST = .5 * EARTH_CIRCUMF

COAST_IS_CLOCKWISE = True

def swap(p):
    return (p[1], p[0])
def antipode(p):
    return (-p[0], p[1]+180.)

def distbear_init(p0):
    v0 = geodesy.ll_to_ecefu(p0)
    vnorth, veast = geodesy.orientate(v0)

    def distbear(p1):
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
    return distbear

def dist_is_axial(dist, cutoff=.1):
    dist = dist % EARTH_FARTHEST
    return (dist < cutoff or dist > EARTH_FARTHEST - cutoff)

class DepthBuffer(object):
    def __init__(self, p0, bearing0, bearing1, res, min_dist):
        self.bearing0 = float(bearing0)
        self.lonspan = (bearing1 - bearing0) % 360. or 360.
        self.size = int(round(self.lonspan / res))
        self.res = self.lonspan / self.size
        self.dist = [-1] * self.size
        self.areas = [set()] * self.size
        self.min_dist = [min_dist] * self.size
        self.max_min_dist = max(self.min_dist)
        self.crossings = [0] * self.size
        self.distbear = distbear_init(p0)
        self.crossings_filled = False

    def output(self):
        self.fill_threshold_crossings()
        return zip(self.dist, self.areas)

    def to_fpx(self, bearing):
        return ((bearing - self.bearing0) % 360.) / self.res

    def to_px(self, bearing):
        px = int(self.to_fpx(bearing))
        assert px >= 0
        return px

    def max_px(self):
        # this is not the size of the viewport, but rather the largest
        # 'virtual' pixel before wraparound. which is the end of the
        # viewport if the range is full 360.
        return self.to_px(self.bearing0 - geodesy.EPSILON)

    def set_dist(self, px, dist, areas):
        if dist < 0:
            dist = self.min_dist[px]
        self.dist[px] = dist
        self.areas[px] = areas

    def post(self, dist, bearing, areas):
        px = self.to_px(bearing)
        if (px < self.size and dist >= self.min_dist[px] and 
            self.dist_threshold(px, dist)):
            self.set_dist(px, dist, areas)
        return px

    def dist_threshold(self, px, dist):
        return self.dist[px] < 0 or dist < self.dist[px]

    def range_test(self, bearing0, bearing1, func):
        fpx0 = self.to_fpx(bearing0)
        fpx1 = self.to_fpx(bearing1)
        px0 = self.to_px(bearing0)
        px1 = self.to_px(bearing1)

        if fpx0 == fpx1:
            # full angle range
            px0, px1 = 0, self.size - 1
        elif fpx0 >= self.size:
            if fpx1 < self.size:
                # overlap-left
                px0 = 0
            elif fpx0 > fpx1:
                # overlap-all
                px0, px1 = 0, self.size - 1
            else:
                # no overlap
                px0, px1 = 0, -1
        elif fpx1 >= self.size:
            # overlap-right
            px1 = self.size - 1
        elif fpx1 < fpx0:
            # wraparound -> overlap-left+right
            px1 = self.size + px1

        return any(func(px % self.size) for px in xrange(px0, px1 + 1))

    def dist_test(self, min_dist, bearing0, bearing1):
        return self.range_test(bearing0, bearing1, lambda px: self.dist_threshold(px, min_dist))

    def threshold_relevancy_test(self, mindist, maxdist, bearing0, bearing1):
        return self.range_test(bearing0, bearing1, lambda px: mindist <= self.min_dist[px] <= maxdist)

    def within_threshold_test(self, maxdist, bearing0, bearing1):
        return not self.range_test(bearing0, bearing1, lambda px: maxdist > self.min_dist[px])
    
    def project_segment(self, p, antip, prev, areas):
        dist, bear, px = self.project_point(p, antip, areas)
        if prev is not None:
            self.fill(p, antip, areas, dist, bear, px, prev)
        return p, dist, bear, px

    def project_point(self, p, antip, areas):
        dist, bear = self.distbear(p)
        if bear is None:
            # exactly antipodal; segment-filling will handle the rest
            bear = 0.
        if antip:
            dist = EARTH_CIRCUMF - dist
            bear = (bear + 180.) % 360.
        px = self.post(dist, bear, areas)
        return dist, bear, px

    def fill(self, p, antip, areas, dist, bear, px, prev):
        # TODO: could bail early if entire segment is out of bearing viewport
        prev_p, prev_dist, prev_bear, prev_px = prev
        adjacent = (abs(px - prev_px) <= 1 or
                    (min(px, prev_px) == 0 and max(px, prev_px) == self.max_px()))
        if adjacent:
            self.check_crossings(dist, px, prev_dist, prev_px, areas)
            return
        if dist_is_axial(dist):
            # segment crosses origin or antipode and will never converge
            # don't bother with crossings as this scenario is rare; we mainly
            # just care about preventing infinite loop
            return

        xyz0 = geodesy.ll_to_ecefu(p)
        xyz1 = geodesy.ll_to_ecefu(prev_p)
        midp = geodesy.ecefu_to_ll(geodesy.vnorm(geodesy.vscale(geodesy.vadd(xyz0, xyz1), .5)))

        middist, midbear, midpx = self.project_point(midp, antip, areas)
        # the ordering of p, mid, and prev is important to maintain direction of coastline
        self.fill(p, antip, areas, dist, bear, px, (midp, middist, midbear, midpx))
        self.fill(midp, antip, areas, middist, midbear, midpx, prev)

    def check_crossings(self, dist, px, prev_dist, prev_px, areas):
        if px >= self.size or prev_px >= self.size:
            # although the crossing logic applies to out-of-view bearings,
            # we don't store a distance threshold, and there is no relevant
            # crossing that could happen at the edge anyway
            return

        def rel(x, px):
            thresh = self.min_dist[px]
            if x > thresh:
                return 1
            elif x < thresh:
                return -1
            else:
                return 0

        rel_a = rel(prev_dist, prev_px)
        rel_b = rel(dist, px)
        if rel_a > rel_b:
            rel_a, rel_b = rel_b, rel_a
            i_close, i_far = 1, 0
        else:
            i_close, i_far = 0, 1

        if (rel_a, rel_b) in ((-1, -1), (1, 1), (0, 1)):
            return
        elif (rel_a, rel_b) == (-1, 1):
            self.add_crossing([prev_px, px][i_close], i_close > i_far, areas)
        elif (rel_a, rel_b) == (-1, 0):
            self.add_crossing([prev_px, px][i_far], i_close > i_far, areas)
        elif (rel_a, rel_b) == (0, 0):
            self.add_crossing(prev_px, True, areas)
            self.add_crossing(px, False, areas)
        else:
            assert False, "can't happen"

    def add_crossing(self, px, is_entry, areas):
        self.set_dist(px, -1, areas)
        # 1 means start fill; -1 means end fill
        self.crossings[px] += (1 if is_entry else -1) * (-1 if COAST_IS_CLOCKWISE else 1)

    def fill_threshold_crossings(self):
        if self.crossings_filled:
            return

        segments = []
        try:
            for i, xing in enumerate(self.crossings):
                if xing == 1:
                    if segments and len(segments[-1]) != 2:
                        raise ValueError
                    segments.append([i])
                elif xing == -1:
                    if segments and len(segments[-1]) != 1:
                        raise ValueError
                    if not segments:
                        segments.append([0, i])
                    else:
                        segments[-1].append(i)
            if segments and len(segments[-1]) == 1:
                segments[-1].append(self.size)
        except ValueError:
            print 'crossings integrity error!'
            for i, c in enumerate(self.crossings):
                if c != 0:
                    print {1: 'start', -1: 'stop'}[c], self.bearing0 + self.res * i, i
            return

        for start, end in segments:
            mid = .5 * (start + end)
            for i in xrange(start, end):
                if start == 0:
                    areas_ix = end
                elif end == self.size:
                    areas_ix = start
                else:
                    areas_ix = start if i < mid else end
                self.set_dist(i, -1, self.areas[areas_ix])

        self.crossings_filled = True

class SegmentSequencer(object):
    def __init__(self, data, p0, db):
        self.data = data
        self.p0 = p0
        self.db = db
        self.q = Queue.PriorityQueue()
        self.all_ixs = set(itertools.chain(*map(pd.ix_ancestor, data.keys())))
        for ix in ('0', '1'):
            for antip in (True, False):
                self.add(ix, antip)

    def add(self, ix, antip):
        bbox = pd.Quad.from_ix(ix)
        bounds = bounds_for_bbox(self.p0, (bbox.y1, bbox.x0), (bbox.y0, bbox.x1))
        if antip:
            bounds = (EARTH_CIRCUMF - bounds[1],
                      EARTH_CIRCUMF - bounds[0],
                      geodesy.anglenorm(bounds[2]+180.),
                      geodesy.anglenorm(bounds[3]+180.))
        self.q.put((bounds[0], ix, antip, bounds))

    def process_next(self):
        _, ix, antip, bounds = self.q.get(False) # empty caught in caller
        mindist, maxdist, minbear, maxbear = bounds

        if mindist > self.db.max_min_dist:
            self.db.fill_threshold_crossings()
        if maxdist < self.db.max_min_dist and self.db.within_threshold_test(maxdist, minbear, maxbear):
            return

        if (not self.db.dist_test(mindist, minbear, maxbear) and
            # include all segments that straddle the threshold, to ensure all crossings are registered
            not self.db.threshold_relevancy_test(mindist, maxdist, minbear, maxbear)):
            return

        print ('%s%s' % (ix, '-' if antip else '+')).ljust(15), mindist
        if ix in self.data:
            return (ix, antip)
        else:
            for i in xrange(4):
                child_ix = ix + str(i)
                if child_ix in self.all_ixs:
                    self.add(child_ix, antip)

    def __iter__(self):
        while True:
            try:
                next = self.process_next()
            except Queue.Empty:
                break
            if next is not None:
                yield next

def postings(p0, data, bearing_res, bearing0, bearing1, near_dist):
    db = DepthBuffer(p0, bearing0, bearing1, bearing_res, near_dist)
    for ix, antip in SegmentSequencer(data, p0, db):
        segment = data[ix]
        for line, areas in segment.iteritems():
            coords = line.coords
            # a ring duplicates the first coord as its last coord, so we don't
            # need to close the loop manually

            prev = None
            for v in coords:
                prev = db.project_segment(swap(v), antip, prev, areas)

    return db.output()




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
        ortho = geodesy.vnorm(geodesy.crossp(*map(geodesy.vnorm, vects)))
    center = geodesy.ecefu_to_ll(ortho)
    dist = geodesy.distance(center, (latmin, lonmin))
    return center, dist

def bounds_to_circle(p0, center, radius):
    dist = geodesy.distance(p0, center)
    mindist = max(0, dist - radius)
    maxdist = min(EARTH_FARTHEST, dist + radius)

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


def lat_for_isodist_tangent(p0, lon):
    dlon = geodesy.anglenorm(lon - p0[1])
    x, y, z = geodesy.ll_to_ecefu((p0[0], dlon))
    if abs(x) < geodesy.EPSILON:
        return None
    if x < 0:
        x = -x
        z = -z
    return math.degrees(math.atan2(z, x))

# inclusive
def in_lon_range(lon, lonmin, lonmax):
    dlon = (lon - lonmin) % 360.
    return dlon >= 0 and dlon <= (lonmax - lonmin) % 360.

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

    dists = [geodesy.distance(p0, key_point) for key_point in key_points()]

    if lat >= latmin and lat <= latmax and in_lon_range(lon, lonmin, lonmax):
        #log('min at pode')
        mindist = 0
    else:
        #log('min at %s' % str(min(key_points(), key=lambda e: geodesy.distance(p0, e))))
        mindist = min(dists)

    if lat >= anti_latmin and lat <= anti_latmax and in_lon_range(lon, anti_lonmin, anti_lonmax):
        #log('max at antipode')
        maxdist = .5 * EARTH_CIRCUMF
    else:
        #log('max at %s' % str(max(key_points(), key=lambda e: geodesy.distance(p0, e))))
        maxdist = max(dists)

    return mindist, maxdist

def bounds_for_bbox(p0, ll, ur):
    center, radius = enclosing_circle_for_bbox(ll, ur)
    mindist, maxdist = minmax_bbox_dist(p0, ll, ur)
    _, _, minbear, maxbear = bounds_to_circle(p0, center, radius)
    return mindist, maxdist, minbear, maxbear
    
