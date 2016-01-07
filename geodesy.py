import math

#geodesic computations, such as distances and bearings between points
#all calculations assume a spherical earth for now

EARTH_EQ_RAD = 6378137.0
EARTH_POL_RAD = 6356752.3
EARTH_MEAN_RAD = 6371009.0

EPSILON = 1.0e-9

#return i dot j
#i and j are any vector of the same dimension
def dotp(i, j):
    return sum([x * y for (x, y) in zip(i, j)])

#return i cross j
#i and j are any 3-vector
def crossp(i, j):
    indexes = [(k, (k + 1) % 3) for k in (1, 2, 0)]
    return tuple([i[t] * j[u] - i[u] * j[t] for (t, u) in indexes])

#return norm/length of v
#v is any vector in any dimension
def vlen(v):
    return math.sqrt(sum([k**2. for k in v]))

#return k * v
#v is any vector in any dimension, k is a scalar factor
def vscale(v, k):
    return tuple([x * k for x in v])

#return i + j
#i and j are any vector of the same dimension
def vadd(i, j):
    return tuple([x + y for (x, y) in zip(i, j)])

#return i - j
#i and j are any vector of the same dimension
def vdiff(i, j):
    return vadd(i, vscale(j, -1.))

#normalize v; exception if norm(v) == 0
#v is any vector in any dimension
def vnorm(v):
    norm = vlen(v)
    if (norm < EPSILON):
        raise ZeroDivisionError()
    return vscale(v, 1. / norm)

#return component of j orthogonal to i, and cosine of angle between i and j
#i and j are unit 3-vectors
def vortho(i, j):
    kcos = dotp(i, j)
    return (vadd(j, vscale(i, -kcos)), kcos)

#rotate v by angle theta around axis, clockwise for positive theta when looking
#from 'axis' toward origin. 
#v and axis are unit 3-vectors, v and axis need not be orthogonal, theta in radians
def vrot(v, axis, theta):
    return vroter(v, axis)(theta)

#vrot, but for many thetas at once
def vroter(v, axis):
    (vo, kcos) = vortho(axis, v)
    vaxial = vscale(axis, kcos)
    vd = crossp(vo, axis)
    return lambda theta: vadd(vaxial, vangle(vo, vd, theta))

#return point on unit sphere corresponding to (lat, lon)
def ll_to_ecefu((lat, lon)):
    rlat = math.radians(lat)
    rlon = math.radians(lon)
    latcos = math.cos(rlat)
    return (math.cos(rlon) * latcos, math.sin(rlon) * latcos, math.sin(rlat))

#convert a point on the unit sphere to (lat, lon)
def ecefu_to_ll((x, y, z)):
    rlat = math.asin(clamp(z, -1., 1.))
    if (abs(x) < EPSILON and abs(y) < EPSILON):
        rlon = 0.
    else:
        rlon = math.atan2(y, x)
    return (math.degrees(rlat), math.degrees(rlon))

#return 'north' and 'east' vectors for a given position vector
def orientate(vp):
    try:
        veast = vnorm(crossp((0., 0., 1.), vp))
    except ZeroDivisionError:
        #at a pole
        veast = (0., -vp[2], 0.)
    vnorth = crossp(vp, veast)
    return (vnorth, veast)

#create angle vector given orthogonal basis vectors 'u' and 'v', and angle 'theta' in radians
def vangle(u, v, theta):
    return vadd(vscale(u, math.cos(theta)), vscale(v, math.sin(theta)))

#return bearing vector for a given position vector and bearing angle
def vbear(vp, bearing):
    rbear = math.radians(bearing)
    (vnorth, veast) = orientate(vp)
    return vangle(vnorth, veast, rbear)

#return distance, in meters, between lat/lon coordinates p0 and p1
def distance(p0, p1, rad=EARTH_MEAN_RAD):
    [v0, v1] = [ll_to_ecefu(p) for p in [p0, p1]]
    (vo, kcos) = vortho(v0, v1)
    ksin = vlen(vo)
    return rad * math.atan2(ksin, kcos)

#return compass bearing from src to dst; None if src/dst are antipodal;
#if src is polar, treat direction of 0 longitude as north
def bearing(src, dst):
    [vsrc, vdst] = [ll_to_ecefu(p) for p in [src, dst]]

    (vdir, _) = vortho(vsrc, vdst)
    if vlen(vdir) < EPSILON:
        #antipodal
        return None

    return _bearing(vsrc, vdir)

def _bearing(vp, vdir):
    (vnorth, veast) = orientate(vp)
    return _xy_to_bearing(dotp(vdir, veast), dotp(vdir, vnorth))

def _xy_to_bearing(x, y):
    return math.degrees(math.atan2(x, y))

#return the coordinates of the position 'distance' meters away from 'p', in direction 'bearing'
#as well as new bearing at the target point 
def plot(p, bearing, distance):
    return plotter_d(p, bearing)(distance)

#plot, but for many distances at once (useful for great circle arcs)
def plotter_d(p, bearing):
    vp = ll_to_ecefu(p)
    vdir = vbear(vp, bearing)
    return lambda d: _plot(vp, vdir, d / EARTH_MEAN_RAD)

#given pos and dir vectors, return plotted position and direction bearing
#after distance 'theta'
def _plot(vp, vdir, theta):
    vp_plot = vangle(vp, vdir, theta)
    vdir_plot = vangle(vdir, vscale(vp, -1.), theta)

    p_plot = ecefu_to_ll(vp_plot)
    bearing_plot = _bearing(vp_plot, vdir_plot)

    return (p_plot, bearing_plot)

#plot, but for many bearings at once (useful for distance range arcs)
def plotter_b(p, distance):
    vp = ll_to_ecefu(p)
    dst = ll_to_ecefu(plot(p, 0, distance)[0])
    rot = vroter(dst, vp)
    return lambda b: ecefu_to_ll(rot(math.radians(b)))


def great_circle(p, bearing, distdelta, distmax=2.*math.pi*EARTH_MEAN_RAD, distmin=0.):
    plot = plotter_v(p, bearing)
    for d in rangef(distmin, distmax, distdelta):
        yield plot(d)[0]

def distance_arc(p, distance, angledelta, anglemin=None, anglemax=None):
    plot = plotter_b(p, distance)
    for a in rangea(angledelta, anglemin, anglemax):
        yield plot(a)


def rangef(min, max, step):
    if feq(step, 0.) and not feq(min, max):
        raise ValueError

    f = min
    while flt(f, max):
        yield f
        f += step
    if feq(min, max) or (flt(min, max) and not feq(f, max + step)):
        yield max

def rangea(step, min=None, max=None):
    normfunc = anglenorm
    if (min == None or max == None or flte(360., max - min)):
        r = rangef(-180., 180., step)
        normfunc = lambda a: a
    elif feq(min, max):
        r = [min]
    elif max > min:
        r = rangef(min, max, step)
    elif flte(min - max, 360.):
        r = rangef(min, max + 360., step)
    else:
        r = []
    for a in r:
        yield normfunc(a)


#offset = 0. - lowest angle
def anglenorm(a, offset=180.):
    return (a + offset) % 360. - offset

def clamp(x, min, max):
    if x > max:
        return max
    elif x < min:
        return min
    else:
        return x

def feq(a, b):
    return abs(a - b) < EPSILON

def flt(a, b):
    return a < b - EPSILON

def flte(a, b):
    return a < b + EPSILON

