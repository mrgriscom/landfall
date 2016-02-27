import sys
import os.path

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

from tornado.ioloop import IOLoop
import tornado.web as web
import tornado.gen as gen
from tornado.template import Template
import tornado.websocket as websocket
from optparse import OptionParser
import logging
import json
import config
import geodesy
import pickle
import landfall as lf
from datetime import datetime
import itertools
import process_data as pd
import random
import math
import util as u

from munsell import munsell as m
m.init()

IX = None
def load_index():
    global IX
    if IX is None:
        print 'loading index...'
        IX = pickle.load(open('data/tmp/tagged_coastline'))
        print 'index loaded.'
        # TODO rebalance index chunks so that each has at least 500 points or so
    return IX

class LandfallHandler(web.RequestHandler):
    def get(self): # technically, should be POST
        _origin = self.get_argument('origin')
        size = int(self.get_argument('size'))
        _range = self.get_argument('range', '0,360')
        mindist = float(self.get_argument('mindist', '100'))

        origin = map(float, _origin.split(',')[:2])
        range = map(float, _range.split(',')[:2])

        assert origin[0] >= -90. and origin[0] <= 90.
        origin[1] = geodesy.anglenorm(origin[1])
        assert size > 0
        range = map(geodesy.anglenorm, range)
        if range[0] == range[1]:
            range[1] += 360.
        assert mindist > 0

        lonspan = (range[1] - range[0]) % 360. or 360.
        res = lonspan / size
        print origin, size, range, lonspan, res, mindist

        print 'generating...'
        postings = lf.project_landfall(origin, load_index(), res, range[0], range[1], mindist)

        print 'saving...'
        output = {
            'origin': origin,
            'res': res,
            'range': range,
            'min_dist': mindist,
            'postings': postings,
        }
        tag = '%s.json' % datetime.now().strftime('%Y%m%d-%H%M%S')
        with open(os.path.join(config.OUTPUT_PATH, tag), 'w') as f:
            json.dump(output, f, indent=2)
        self.redirect(self.reverse_url('render', tag))
        
class OutputListHandler(web.RequestHandler):
    def get(self):
        outputs = []
        for filename in sorted(os.listdir(config.OUTPUT_PATH)):
            with open(os.path.join(config.OUTPUT_PATH, filename)) as f:
                data = json.load(f)
            del data['postings']
            data['file'] = filename
            data['lonspan'] = (data['range'][1] - data['range'][0]) % 360. or 360.
            data['width'] = int(round(data['lonspan'] / data['res']))
            outputs.append(data)
        self.render('list.html', outputs=outputs)

# TODO panorama exif

class RenderHandler(web.RequestHandler):
    def get(self, tag):
        width = int(self.get_argument('width', '3000'))
        height = int(self.get_argument('height', '800'))
        min_downscale = float(self.get_argument('mindownscale', '2'))

        min_dist = float(self.get_argument('mindist', '0'))
        max_dist = float(self.get_argument('maxdist', '0'))

        # custom ticks for distance axis -- stopgap for better automatic logic
        yticks = self.get_argument('yticks', None)
        if yticks:
            yticks = [float(y) if y != 'antip' else y for y in yticks.split(',')]

        num_colors = int(self.get_argument('numcolors', '6'))
        hues = self.get_argument('hues', '')
        if hues:
            hues = map(float, hues.split(','))
        else:
            hues = [int(360. * float(i) / num_colors) for i in xrange(num_colors)]
        lum_close = float(self.get_argument('lumclose', '.5'))
        lum_far = float(self.get_argument('lumfar', '.7'))
        sat_close = float(self.get_argument('satclose', '.6'))
        sat_far = float(self.get_argument('satfar', '.03'))
        def get_ramp(hue):
            key = (hue, lum_close, lum_far, sat_close, sat_far)
            if key not in COLOR_CACHE:
                print 'generating color ramp: hue %s lum %s %s sat %s %s' % key
                COLOR_CACHE[key] = color_ramp(*key)
            return COLOR_CACHE[key]
        colors = map(get_ramp, hues)
        color_near_dist = float(self.get_argument('colorneardist', '0'))
        color_far_dist = float(self.get_argument('colorfardist', '0'))
        assert all(k < lf.EARTH_CIRCUMF for k in (color_near_dist, color_far_dist))

        force_color = self.get_argument('forcecolor', '')
        force_color = dict(e.split(':') for e in force_color.split(',') if e)
        for k in force_color.keys():
            force_color[k] = int(force_color[k])
        force_interfere = self.get_argument('diffcolor', '')
        force_interfere = list(itertools.chain(*(itertools.combinations(pair.split(':'), 2) for pair in force_interfere.split(',') if pair)))
        no_subdivisions = self.get_argument('nosubdiv', '')
        no_subdivisions = set(filter(None, no_subdivisions.split(',')))
        resolve_as = self.get_argument('resolve', '')
        resolve_as = dict(pair.split(':') for pair in resolve_as.split(',') if pair)
        random_seed = self.get_argument('randseed', '')
        random_seed = int(random_seed) if random_seed else None

        dist_unit = self.get_argument('distunit', 'km')
        assert dist_unit in ('km', 'mi', 'nmi', 'deg', 'none')
        if dist_unit == 'none':
            dist_unit = None

        # TODO nosubdiv should modify explicitly specified cc's in params
        params = {
            'dim': [width, height],
            'min_dist': min_dist,
            'max_dist': max_dist,
            'colors': colors,
            'colordists': [color_near_dist, color_far_dist],
            'force_color': force_color,
            'force_interfere': force_interfere,
            'no_subdivisions': list(no_subdivisions),
            'resolve_as': resolve_as,
            'random_seed': random_seed,
            'dist_unit': dist_unit,
            'yticks': yticks,
        }

        with open(os.path.join(config.OUTPUT_PATH, tag)) as f:
            data = json.load(f)
        data['lonspan'] = (data['range'][1] - data['range'][0]) % 360. or 360.
        data['wraparound'] = (data['lonspan'] == 360.)

        left = self.get_argument('left', '')
        left = float(left) if left else data['range'][0]
        right = self.get_argument('right', '')
        right = float(right) if right else data['range'][1]
        trimleft = float(self.get_argument('trimleft', '0'))
        trimright = float(self.get_argument('trimright', '0'))
        if trimleft != 0 or trimright != 0:
            assert not data['wraparound'], 'trim params not valid for full-360 wraparound; consider bearmargin'
        left += trimleft
        right -= trimright
        params['bear0'] = left
        params['bearspan'] = (right - left) % 360. or 360.
        if not data['wraparound']:
            leftrel = (left - data['range'][0]) % 360.
            rightrel = (right - data['range'][0]) % 360.
            assert leftrel < data['lonspan'] and rightrel <= data['lonspan'] and leftrel < rightrel, 'specified range out of data range'
        bearcenter = self.get_argument('bearcenter', '')
        bearmargin = self.get_argument('bearmargin', '')
        if bearcenter or bearmargin:
            assert data['wraparound'], 'only valid for full-360 wraparound'
            bearcenter = float(bearcenter) if bearcenter else .5*(data['range'][0] + data['range'][1])
            bearmargin = float(bearmargin) if bearmargin else 0
            params['bear0'] = bearcenter - 180. - bearmargin
            params['bearspan'] = 360. + 2. * bearmargin
        params['res'] = float(params['bearspan']) / width

        self.process_downsampling(data, params['res'], min_downscale)
        self.process_admins(data, params)

        self.render('render.html', data=data, params=params)

    def process_downsampling(self, data, bearres, min_downscale):
        data['size'] = len(data['postings'])
        downscale = math.log(bearres / data['res'], 2)
        if downscale < 0:
            raise ValueError('rendering width greater than data size')
        elif downscale < min_downscale:
            print 'WARNING: Max recommended width is %dpx (downscaling factor: %s); suggest recomputing the source data at higher resolution to avoid aliasing' % (math.floor(data['size'] * 2**-min_downscale), min_downscale)
        downsample = 2**max(math.floor(downscale - min_downscale), 0)
        if downsample > 1:
            print 'Downsampling by %dx' % downsample
            new_size = int(round(float(data['size']) / downsample))
            scale_factor = float(new_size) / data['size']
            downsampled = [(-1, None)] * new_size
            for i, posting in enumerate(data['postings']):
                dist, admin = posting
                if dist < 0:
                    continue
                i_ds = int(math.floor(i * scale_factor))
                ds_dist = downsampled[i_ds][0]
                if ds_dist < 0 or dist < ds_dist:
                    downsampled[i_ds] = posting
            data['postings'] = downsampled
            data['size'] = new_size
            data['res'] *= downsample

    def process_admins(self, data, params):
        from pprint import pprint #debug

        admin_postings = [(p[1] or pd.CC_UNCLAIMED) for p in data['postings']]

        def _admin_info():
            # load all admin info
            info = {}
            admin_info = pd.load_admin_info(pd.admin_info_path)
            for k, v in admin_info.iteritems():
                info[k] = {
                    'name': v['name'] + (' (%s)' % v['parent'] if v['parent'] else ''),
                    'parent': v['parent'] if v['type'] == 'subdivision' else None,
                }
            for k, v in config.disputed_areas.iteritems():
                info['%s-%s' % (pd.CC_DISPUTED, k)] = {'name': v['name'], 'parent': None}
            info[pd.CC_UNCLAIMED] = {'name': 'terra nullius', 'parent': None}

            # filter to relevant admins:
            # admins seen in output
            entities = set(admin_postings)
            # admins that could be converted to
            entities.update(params['resolve_as'].values())
            # and parents thereof
            parents = set(filter(None, (info[e]['parent'] for e in entities)))
            entities.update(parents)
            return dict((k, v) for k, v in info.iteritems() if k in entities)
        admin_info = _admin_info()
        data['admin_info'] = admin_info

        def change_admin(a):
            # resolution may be a subdivision so handle before checking subdiv suppression
            a = params['resolve_as'].get(a, a)
            parent = admin_info[a]['parent']
            if parent and parent in params['no_subdivisions']:
                a = parent
            return a
        admin_postings = map(change_admin, admin_postings)
        data['admin_postings'] = admin_postings
        admins = set(admin_postings)
        print '\n'.join(sorted('%s %s' % (a, admin_info[a]['name']) for a in admins))

        def ix_offset(i, offset, size, wrap=None):
            # wrap=None (default) = use wrap value of data
            # wrap=T/F = force wrap mode
            wrap = data['wraparound'] if wrap is None else wrap
            j = i + offset
            if wrap:
                return j % size
            else:
                return j if 0 <= j < size else -1

        def _admin_segments():
            segments = [];
            segment = None
            for i in xrange(data['size']):
                inext = ix_offset(i, 1, data['size'])
                if inext < 0:
                    break

                admin0 = admin_postings[i]
                admin1 = admin_postings[inext]
                if i == 0:
                    segment = {'admin': admin0, 'start': 0}
                if admin0 != admin1:
                    segment['end'] = i
                    segments.append(segment)
                    segment = {'admin': admin1, 'start': inext}
            if data['wraparound']:
                if segments:
                    segments[0]['start'] = segment['start']
                else:
                    # no transitions
                    pass
            else:
                segment['end'] = data['size'] - 1
                segments.append(segment)
            return segments
        admin_segments = _admin_segments()

        interferences = {}
        interf_prio = []
        def interfere(a, b, type=None):
            type = type or interf_prio[-1]
            if a != b:
                key = tuple(sorted((a, b)))
                if key not in interferences:
                    interferences[key] = type
        
        interf_prio.append('force')
        for a, b in params['force_interfere']:
            assert a in admins, a
            assert b in admins, b
            interfere(a, b)

        for cc, i in params['force_color'].iteritems():
            assert cc in admins, cc
            assert 0 <= i < len(params['colors']), i

        def search_adjacent(i, forward, force_wrap=False):
            dir = 1 if forward else -1
            for delta in xrange(1, len(admin_segments)):
                i_adj = ix_offset(i, dir*delta, len(admin_segments), True if force_wrap else None)
                if i_adj < 0:
                    break
                yield admin_segments[i_adj]

        MAX_ADJ_INTERF = 2  # TODO make param (necessary?)
        for delta in xrange(1, MAX_ADJ_INTERF + 1):
            interf_prio.append('adj%d' % delta)
            for i in xrange(len(admin_segments)):
                admin = admin_segments[i]['admin']
                for dir in (True, False):
                    seen = set([admin])
                    for adj in search_adjacent(i, dir):
                        adj_admin = adj['admin']
                        seen.add(adj_admin)
                        if len(seen) == delta + 1:
                            interfere(admin, adj_admin)
                            break

        MAX_PX_INTERF = 20  # screen pixels; TODO make param
        MAX_POSTINGS_INTERF = int(round(float(MAX_PX_INTERF) * params['res'] / data['res']))
        interf_prio.append('nearpx')
        world_width = int(round(360. / data['res']))
        assert world_width >= data['size']
        assert not data['wraparound'] or world_width == data['size']
        for i, seg in enumerate(admin_segments):
            for forward in (True, False):
                for adj in search_adjacent(i, forward, force_wrap=True):
                    if forward:
                        if (adj['start'] - seg['end']) % world_width > MAX_POSTINGS_INTERF:
                            break
                    else:
                        if (seg['start'] - adj['end']) % world_width > MAX_POSTINGS_INTERF:
                            break
                    interfere(seg['admin'], adj['admin'])

        #pprint(interferences)
        #print len(interferences)

        costs = dict((type, 50**(len(interf_prio) - 1 - i)) for i, type in enumerate(interf_prio))
        adjacency = dict((edge, costs[type]) for edge, type in interferences.iteritems())
        adjacency.update((tuple(reversed(edge)), costs[type]) for edge, type in interferences.iteritems())
        adjacent_to = u.map_reduce(adjacency.keys(), lambda edge: [edge], set)

        # INITIALIZE RANDOMNESS
        seed = params['random_seed']
        if seed is None:
            seed = hash(random.getstate())
        random.seed(seed)
        print 'random seed: %s' % seed

        def rand_color():
            return random.randint(0, len(params['colors']) - 1)
        # initial state is random
        colors = dict((a, params['force_color'].get(a, rand_color())) for a in admins)
        color_keys = list(set(colors.keys()) - set(params['force_color']))
        def coloring_energy(colors):
            return sum(cost for edge, cost in adjacency.iteritems() if edge[0] < edge[1] and colors[edge[0]] == colors[edge[1]])
        def energy_diff(colors, move):
            def node_diff_cost(node, newcolor, ignore=[]):
                def node_cost(color):
                    return sum(adjacency[(node, neighbor)] for neighbor in adjacent_to[node] if neighbor not in ignore and color == colors[neighbor])
                return node_cost(newcolor) - node_cost(colors[node])

            if move[0] == 'change':
                cc, color = move[1:]
                if colors[cc] == color:
                    return 0
                return node_diff_cost(cc, color)
            elif move[0] == 'swap':
                a, b = move[1:]
                if colors[a] == colors[b]:
                    return 0
                return node_diff_cost(a, colors[b], [b]) + node_diff_cost(b, colors[a], [a])
        def gen_move():
            if random.random() < .5:
                # change color
                return ('change', random.choice(color_keys), rand_color())
            else:
                # swap colors
                pair = random.sample(color_keys, 2) if len(color_keys) >= 2 else (0, 0)
                return ('swap', pair[0], pair[1])
        def apply_move(colors, move):
            if move[0] == 'change':
                cc, color = move[1:]
                colors[cc] = color
            elif move[0] == 'swap':
                a, b = move[1:]
                colors[a], colors[b] = colors[b], colors[a]
            return colors

        # baseline speed of cooling -- calibrated so that each node is touched ~10 times per 10% drop in temperature
        COOLING_BASELINE = 0.99

        # average # of nodes touched per move (50% change color (1 node) + 50% swap (2 nodes))
        AVG_NODES_TOUCHED_PER_MOVE = 1.5

        # failsafe to stop searching -- should normally terminate earlier by reaching a frozen state
        MIN_TEMPERATURE = 0.0001

        # cooling amount per iteration
        cooling_factor = COOLING_BASELINE**(AVG_NODES_TOUCHED_PER_MOVE / len(color_keys))

        # if no change to energy (within 'frozen_window') after this much relative temperature drop, terminate
        frozen_threshold = .9
        frozen_window = .5 # since all our costs are integers, just has to be less than 1

        # set initial temperature where this proportion of worse moves are accepted (conventional wisdom suggests
        # .5 but this problem space seems to converge well from a lower starting temperature)
        init_worse_accept_p = .1

        current_energy = coloring_energy(colors)

        temp_baseline_iterations = 100
        temp_baseline_failsafe_iterations = 1000
        baselines = []
        for i in xrange(temp_baseline_failsafe_iterations):
            ediff = energy_diff(colors, gen_move())
            if ediff > 0:
                baselines.append(ediff)
                if len(baselines) == temp_baseline_iterations:
                    break

        def intercept(fn, min, max, res):
            mid = .5*(min + max)
            if max - min < res:
                return mid
            elif fn(mid) < 0:
                return intercept(fn, mid, max, res)
            else:
                return intercept(fn, min, mid, res)
        def accept_rate(temp):
            if baselines:
                return sum(math.exp(-ediff / temp) for ediff in baselines) / len(baselines)
            else:
                return 0. # just a failsafe in case we can't get any worse examples to calibrate
                          # which would imply the current state is pretty fucking bad, thus force
                          # the temperature to the max possible. alternatively, the state space
                          # may be very small, in which case the 'frozen' state will occur
                          # immediately.
        temperature = math.exp(intercept(lambda x: accept_rate(math.exp(x)) - init_worse_accept_p, 0., math.log(1e9), .1))
        print 'init temperature:', temperature

        i = 0
        frozen_at = temperature
        frozen_energy = current_energy
        while temperature > MIN_TEMPERATURE:
            if current_energy == 0:
                break
            if temperature / frozen_at < frozen_threshold:
                break

            move = gen_move()
            ediff = energy_diff(colors, move)

            if ediff <= 0.:
                acceptance_p = 1.
                accept = True
                state = 'better'
            else:
                acceptance_p = math.exp(-ediff / temperature)
                accept = (random.random() < acceptance_p)
                state = 'worse' if accept else 'ign'
            #print '%.5f % 8d % 8d %.5f %s' % (temperature, current_energy, current_energy+ediff, acceptance_p, state)

            if accept:
                colors = apply_move(colors, move)
                current_energy += ediff
                if abs(frozen_energy - current_energy) > frozen_window:
                    frozen_at = temperature
                    frozen_energy = current_energy

            i += 1
            temperature *= cooling_factor
        print '%d simulated annealing iterations' % i

        data['colors'] = colors
        #print data['colors']

        print 'conflicts:'
        for edge, type in interferences.iteritems():
            if colors[edge[0]] == colors[edge[1]]:
                print edge, type




class KmlHandler(web.RequestHandler):
    def get(self, tag):
        return None

        """
        with open(os.path.join(config.OUTPUT_PATH, tag)) as f:
            data = json.load(f)

        def get_bearing(i):
            return data['range'][0] + float(i) * data['res']

        def contiguous_segments():
            start = 0
            for i in xrange(1, data['size']):
                j = (i + 1) % data['size']
                if j == 0 and not data['wraparound']:
                    break
                dist_i = data['postings'][i]
                dist_j = data['postings'][j]
                quantum = min(dist_i, dist_j) * math.radians(data['res'])
                if abs(dist_i - dist_j) >= 10. * quantum:
                    yield (start, i)
                    start = j
            if data['wraparound']:
                # shit, already yielded
            else:
                yield (start, data['size'] - 1)


        # wraparound or constrained
        # split into contiguous segments (edges vs. centers)
        # split by admin
        # max # points

        # discontinuities
        # max ray length

        # mindist and bear limits

        # TODO(clean this up):

        antipode = (-data['origin'][0], geodesy.anglenorm(data['origin'][1] + 180.))

        path = []
        prevdist = None
        for i, (dist, _) in enumerate(itertools.chain(data['postings'], [data['postings'][0]])):
            if dist < 0:
                dist = lf.EARTH_CIRCUMF - data['min_dist'] # TODO go all the way back to start but insert interstitial point
            bear = data['range'][0] + i * data['res']
            p = geodesy.plot(data['origin'], bear, dist)[0]
            if prevdist is not None and (dist > lf.EARTH_FARTHEST) != (prevdist > lf.EARTH_FARTHEST):
                path.append(antipode)
            path.append(p)
            prevdist = dist

        SEGSIZE = 1000
        segments = [path[i:min(i+1+SEGSIZE, len(path))] for i in xrange(0, len(path), SEGSIZE)]
        segments.append([geodesy.plot(data['origin'], bear, data['min_dist'])[0] for bear in xrange(0, 361, 5)])

        self.set_header('Content-Type', 'application/vnd.google-earth.kml+xml')
        self.set_header('Content-Disposition', 'attachment; filename="landfall.kml"')
        self.render('render.kml', segments=segments)
"""

COLOR_CACHE = {}

def color_ramp(hue, lummin, lummax, satmin, satmax, num_steps=1000):
    hue /= 360.

    def get_color(hue, lum, sat):
        rgb = m.munsell(hue, lum, sat)
        if not m.in_gamut(rgb):
            # binary search to find most saturated in-gamut color for hue/lum
            sat_quantum = .001
            sat_lo = 0.
            sat_hi = sat
            while sat_hi - sat_lo > sat_quantum:
                sat = .5 * (sat_lo + sat_hi)
                cand = m.munsell(hue, lum, sat)
                if m.in_gamut(cand):
                    rgb = cand
                    sat_lo = sat
                else:
                    sat_hi = sat
        return '#%s' % ''.join('%02x' % k for k in m.rgb_to_hex(rgb))

    def color_for_k(k):
        lum = lummin * (1-k) + lummax * k
        sat = satmin * (1-k) + satmax * k
        return get_color(hue, lum, sat)

    steps = [float(i) / (num_steps - 1) for i in xrange(num_steps)]
    return map(color_for_k, steps)

if __name__ == "__main__":

    if not os.path.exists(config.OUTPUT_PATH):
        os.mkdir(config.OUTPUT_PATH)

    try:
        port = int(sys.argv[1])
    except IndexError:
        port = 8000

    application = web.Application([
        (r'/', OutputListHandler),
        (r'/render/(?P<tag>.*)', RenderHandler, {}, 'render'),
        (r'/kml/(?P<tag>.*)', KmlHandler, {}, 'kml'),
        (r'/landfall', LandfallHandler),
        (r'/(.*)', web.StaticFileHandler, {'path': 'static'}),
    ], template_path='templates', debug=True)
    application.listen(port)

    try:
        IOLoop.instance().start()
    except KeyboardInterrupt:
        pass
    except Exception, e:
        print e
        raise

    logging.info('shutting down...')
