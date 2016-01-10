import math
import geodesy
from shapely.geometry import LineString, LinearRing
import process_data as pd

EARTH_CIRCUMF = 2*math.pi*geodesy.EARTH_MEAN_RAD
EARTH_FARTHEST = .5 * EARTH_CIRCUMF

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

class DepthBuffer(object):
    def __init__(self, p0, bearing0, bearing1, res, min_dist):
        self.bearing0 = float(bearing0)
        self.lonspan = (bearing1 - bearing0) % 360. or 360.
        self.size = int(round(self.lonspan / res))
        self.res = self.lonspan / self.size
        self.dist = [-1] * self.size
        self.areas = [set()] * self.size
        self.min_dist = [min_dist] * self.size
        self.distbear = distbear_init(p0)

    def output(self):
        return zip(self.dist, self.areas)

    def to_fpx(self, bearing):
        return ((bearing - self.bearing0) % 360.) / self.res

    def to_px(self, bearing):
        px = int(self.to_fpx(bearing))
        assert px >= 0
        return px

    def post(self, dist, bearing, areas):
        px = self.to_px(bearing)
        if (px < self.size and dist >= self.min_dist[px] and 
            self.dist_threshold(px, dist)):
            self.dist[px] = dist
            self.areas[px] = areas
        return px

    def dist_threshold(self, px, dist):
        return self.dist[px] < 0 or dist < self.dist[px]

    def dist_test(self, min_dist, bearing0, bearing1):
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

        return any(self.dist_threshold(px % self.size, min_dist) for px in xrange(px0, px1 + 1))

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

    def fill(self, p, antip, areas, dist, bear, px, prev, depth=0):
        prev_p, prev_dist, prev_bear, prev_px = prev
        adjacent = (abs(px - prev_px) <= 1 or
                    (min(px, prev_px) == 0 and max(px, prev_px) == self.size - 1 and self.lonspan == 360.))
        if adjacent:
            return
        if depth == 30:
            # lines that cross origin or antipode will never converge
            # TODO find a more elegant way of handling this, likely at same time as handling
            # min-dist crossovers
            return

        xyz0 = geodesy.ll_to_ecefu(p)
        xyz1 = geodesy.ll_to_ecefu(prev_p)
        midp = geodesy.ecefu_to_ll(geodesy.vnorm(geodesy.vscale(geodesy.vadd(xyz0, xyz1), .5)))

        middist, midbear, midpx = self.project_point(midp, antip, areas)
        self.fill(midp, antip, areas, middist, midbear, midpx, (p, dist, bear, px), depth=depth+1)
        self.fill(midp, antip, areas, middist, midbear, midpx, prev, depth=depth+1)

def postings(p0, data, bearing_res, bearing0, bearing1, near_dist):
    db = DepthBuffer(p0, bearing0, bearing1, bearing_res, near_dist)

    thresholds = {}
    for ix in data.keys():
        bbox = pd.Quad.from_ix(ix)
        bounds = bounds_for_bbox(p0, (bbox.y1, bbox.x0), (bbox.y0, bbox.x1))
        thresholds[(ix, False)] = bounds
        thresholds[(ix, True)] = (EARTH_CIRCUMF - bounds[1],
                                  EARTH_CIRCUMF - bounds[0],
                                  geodesy.anglenorm(bounds[2]+180.),
                                  geodesy.anglenorm(bounds[3]+180.))
    segment_order = sorted(thresholds.keys(), key=lambda ix: thresholds[ix][0])

    for segnum, (ix, antip) in enumerate(segment_order):
        mindist, _, minbear, maxbear = thresholds[(ix, antip)]
        if not db.dist_test(mindist, minbear, maxbear):
            continue
        print ix, mindist

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

def bounds_for_bbox(p0, ll, ur):
    center, radius = enclosing_circle_for_bbox(ll, ur)
    # TODO: use exact computation for distance thresholds
    return bounds_to_circle(p0, center, radius)
