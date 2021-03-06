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
from functools import partial

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
    DEFAULT_ANNOTATIONS = {
        'osm': u'data \xa9 OpenStreetMap contributors',
        'mrgris': 'mrgris.com/projects/landfall',
    }

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
        ylabelreps = int(self.get_argument('ylabelrepeat', '1'))
        annot1 = self.get_argument('attrib1', self.DEFAULT_ANNOTATIONS['osm'])
        annot2 = self.get_argument('attrib2', self.DEFAULT_ANNOTATIONS['mrgris'])

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
            'ylabelreps': ylabelreps,
            'primary_annotation': annot1,
            'secondary_annotation': annot2,
        }

        with open(os.path.join(config.OUTPUT_PATH, tag)) as f:
            data = json.load(f)
        data['tag'] = os.path.splitext(tag)[0]
        data['lonspan'] = (data['range'][1] - data['range'][0]) % 360. or 360.
        data['wraparound'] = (data['lonspan'] == 360.)
        data['size'] = len(data['postings'])

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

        # send to client as string or else js will lose precision
        params['random_seed'] = str(params['random_seed'])
        self.render('render.html', data=data, params=params)

    def process_downsampling(self, data, bearres, min_downscale):
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
        params['random_seed'] = seed
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
                    return sum(adjacency[(node, neighbor)] for neighbor in adjacent_to.get(node, []) if neighbor not in ignore and color == colors[neighbor])
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
                try:
                    pair = random.sample(color_keys, 2)
                except ValueError:
                    pair = (color_keys[0], color_keys[0])
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

        # if no change to energy (within 'frozen_window') after this much relative temperature drop, terminate
        frozen_threshold = .9
        frozen_window = .5 # since all our costs are integers, just has to be less than 1

        current_energy = coloring_energy(colors)
        if not color_keys:
            # all colors assigned; abort immediately
            temperature = 0
        else:
            # cooling amount per iteration
            cooling_factor = COOLING_BASELINE**(AVG_NODES_TOUCHED_PER_MOVE / len(color_keys))

            # set initial temperature where this proportion of worse moves are accepted (conventional wisdom suggests
            # .5 but this problem space seems to converge well from a lower starting temperature)
            init_worse_accept_p = .1

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
        with open(os.path.join(config.OUTPUT_PATH, tag)) as f:
            data = json.load(f)
        segments = vector_segments(data)

        def color_for_style(style):
            BOUNDS_COLOR = '8866ff'
            LANDFALL_COLOR = 'ff0066'
            TANGENT_COLOR = 'aaaaff'
            return {
                'bounds': BOUNDS_COLOR,
                'tangent': TANGENT_COLOR,
                'landfall': LANDFALL_COLOR,
            }[style]

        def set_color(seg):
            color = color_for_style(seg['style'])
            if len(color) == 6:
                color += 'ff'
            kmlcolor = ''.join(reversed([color[2*k:2*(k+1)] for k in xrange(4)]))
            return {'color': kmlcolor, 'postings': seg['postings']}
        segments = map(set_color, segments)

        self.set_header('Content-Type', 'application/vnd.google-earth.kml+xml')
        self.set_header('Content-Disposition', 'attachment; filename="landfall.kml"')
        self.render('render.kml', segments=segments, origin=data['origin'])
        
class GeojsonHandler(web.RequestHandler):
    def get(self, tag):
        with open(os.path.join(config.OUTPUT_PATH, tag)) as f:
            data = json.load(f)
        origin = tuple(reversed(data['origin']))
        segments = vector_segments(data)

        fragments = []
        for seg in segments:
            for a, b in pairwise_postings(seg['postings']):
                lon0 = min(a[1], b[1])
                lon1 = max(a[1], b[1])
                if abs(lon1 - lon0) < 180:
                    fragments.append((lon0, lon1))
                else:
                    fragments.append((lon1, 180))
                    fragments.append((-180, lon0))
        fragments.sort()
        gaps = []
        frontier = -180
        for f in fragments:
            if f[0] > frontier:
                gaps.append((frontier, f[0]))
            frontier = max(frontier, f[1])
        if frontier < 180:
            if gaps[0][0] == -180:
                gaps[0] = (frontier - 360, gaps[0][1])
            else:
                gaps.append(frontier, 180)
        if not gaps:
            duplicate = True
        else:
            duplicate = False
            biggest_gap = max(gaps, key=lambda e: e[1] - e[0])
            center = geodesy.anglenorm(.5*sum(biggest_gap) + 180)
        
        segments = transform_segments(segments, partial(make_mercator_safe, tolerance=100))
        
        if duplicate:
            origins = [(origin[0] + 360*i, origin[1]) for i in xrange(-1, 2)]
            duplicated_segments = split_segments(segments, partial(clip_to_window, lon_center=data['origin'][1] - 180))
            duplicated_segments.extend(split_segments(segments, partial(clip_to_window, lon_center=data['origin'][1] + 180)))
            segments = duplicated_segments
        else:
            origins = [(geodesy.anglenorm(origin[0], 180 - center), origin[1])]
            segments = split_segments(segments, partial(clip_to_window, lon_center=center))
            
        def color_for_style(style):
            BOUNDS_COLOR = '8866ff'
            LANDFALL_COLOR = 'ff0066'
            TANGENT_COLOR = 'ffffaa'
            return {
                'bounds': (BOUNDS_COLOR, 1.),
                'tangent': (TANGENT_COLOR, .5),
                'landfall': (LANDFALL_COLOR, 1.),
            }[style]

        def to_feature(segment):
            coords = [(lon, lat) for lat, lon in segment['postings']]
            color, opacity = color_for_style(segment['style'])
            return {
                "type": "Feature",
                "geometry": {
                    "type": "LineString",
                    "coordinates": coords,
                },
                "properties": {
                    "stroke": '#%s' % color,
                    "stroke-opacity": opacity,
                    "stroke-width": 1,
                },
            }

        def make_origin(o):
            return {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": o,
                },
                "properties": {
                    "title": "vantage point",
                    "produced_by": "http://mrgris.com/projects/landfall",
                    "attribution": u'data \xa9 OpenStreetMap contributors',
                },
            }
        
        geojson = {
            "type": "FeatureCollection",
            "features": map(make_origin, origins) + map(to_feature, segments),
        }

        self.set_header('Access-Control-Allow-Origin', '*')        
        self.set_header('Content-Disposition', 'attachment; filename="landfall.geojson"')
        self.write(geojson)
        
def vector_segments(data):
    data['lonspan'] = (data['range'][1] - data['range'][0]) % 360. or 360.
    data['wraparound'] = (data['lonspan'] == 360.)
    data['size'] = len(data['postings'])

    def get_bearing(i):
        return data['range'][0] + float(i) * data['res']

    def ix_offset(i, offset, size, wrap=None):
        # wrap=None (default) = use wrap value of data
        # wrap=T/F = force wrap mode
        wrap = data['wraparound'] if wrap is None else wrap
        j = i + offset
        if wrap:
            return j % size
        else:
            return j if 0 <= j < size else -1

    def quantum(dist):
        effective_radius = geodesy.EARTH_MEAN_RAD * abs(math.sin(dist / geodesy.EARTH_MEAN_RAD))
        return effective_radius * math.radians(data['res'])

    def _contig_segments():
        DISCONT_THRESHOLD = 1. / math.tan(math.radians(5.))

        segments = [];
        segment = None
        for i in xrange(data['size']):
            inext = ix_offset(i, 1, data['size'])
            if inext < 0:
                break

            dist0 = data['postings'][i][0]
            dist1 = data['postings'][inext][0]
            hasdist0 = dist0 > 0
            hasdist1 = dist1 > 0
            discontinuity = False
            if hasdist0 != hasdist1:
                discontinuity = True
            elif hasdist0 and hasdist1:
                if abs(dist0 - dist1) > DISCONT_THRESHOLD * min(quantum(dist0), quantum(dist1)):
                    discontinuity = True

            if i == 0:
                segment = {'start': 0}
            if discontinuity:
                segment['end'] = i
                segments.append(segment)
                segment = {'start': inext}
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
    contig_segments = _contig_segments()

    draw_segments = []

    near_dist = []
    for a in geodesy.rangea(3., data['range'][0], data['range'][1]):
        near_dist.append((a, data['min_dist']))
    draw_segments.append({'style': 'bounds', 'postings': near_dist})
    if not data['wraparound']:
        for k in (0, data['size']):
            dist = data['postings'][k if k == 0 else -1][0]
            if dist < 0:
                dist = 2*math.pi*geodesy.EARTH_MEAN_RAD
            bearing = get_bearing(k)
            draw_segments.append({'style': 'bounds', 'postings': [(bearing, data['min_dist']), (bearing, dist)]})

    for i in xrange(len(contig_segments)):
        seg = contig_segments[i]
        if i == len(contig_segments) - 1 and not data['wraparound']:
            break
        nextseg = contig_segments[(i+1) % len(contig_segments)]
        
        dist0 = data['postings'][seg['end']][0]
        dist1 = data['postings'][nextseg['start']][0]
        if dist0 < 0:
            dist0 = dist1 + 2*math.pi*geodesy.EARTH_MEAN_RAD
        elif dist1 < 0:
            dist1 = dist0 + 2*math.pi*geodesy.EARTH_MEAN_RAD
        
        bearing = get_bearing(nextseg['start'])
        draw_segments.append({'style': 'tangent', 'postings': [(bearing, dist0), (bearing, dist1)]})

    contig_segments = filter(lambda seg: data['postings'][seg['start']][0] > 0, contig_segments)
    for seg in contig_segments:
        postings = []
        def wrap_ix():
            for i in xrange((seg['end'] - seg['start']) % data['size'] + 1):
                yield (seg['start'] + i) % data['size']
        for ix in wrap_ix():
            if ix == seg['start']:
                bi = ix
            elif ix == seg['end']:
                bi = ix + 1
            else:
                bi = ix + .5
            postings.append((get_bearing(bi), data['postings'][ix][0]))
        if seg['start'] == seg['end']:
            postings.append((get_bearing(seg['end'] + 1), data['postings'][seg['end']][0]))
        draw_segments.append({'style': 'landfall', 'postings': postings})

    def max_point_spacing(seg_postings):
        MAX_SPACING = .5*math.pi*geodesy.EARTH_MEAN_RAD
        yield seg_postings[0]
        for a, b in pairwise_postings(seg_postings):
            if a[0] != b[0]:
                yield b
                continue
            for k in list(geodesy.rangef(0, abs(a[1] - b[1]), MAX_SPACING))[1:]:
                yield (b[0], a[1] + (1 if b[1] > a[1] else -1) * k)
    draw_segments = transform_segments(draw_segments, max_point_spacing)

    def project_segment(postings):
        return [geodesy.plot(data['origin'], bear, dist)[0] for (bear, dist) in postings]
    draw_segments = transform_segments(draw_segments, project_segment)

    def max_segment_points(postings):
        MAX_POINTS_PER_SEGMENT = 2000
        num_subsegs = int(math.ceil((len(postings) - 1.) / (MAX_POINTS_PER_SEGMENT - 1.)))
        for i in xrange(num_subsegs):
            start = i * (MAX_POINTS_PER_SEGMENT - 1)
            yield postings[start:start+MAX_POINTS_PER_SEGMENT]
    draw_segments = split_segments(draw_segments, max_segment_points)

    return draw_segments

def pairwise_postings(postings):
    for i in xrange(len(postings) - 1):
        yield (postings[i], postings[i + 1])

def copy_segment(segment, new_postings):
    new_seg = dict(segment)
    new_seg['postings'] = new_postings
    return new_seg

def transform_segments(segments, transform):
    return map(lambda seg: copy_segment(seg, list(transform(seg['postings']))), segments)

def split_segments(segments, split):
    def split_segment(segment):
        for postings in split(segment['postings']):
            yield copy_segment(segment, postings)
    return list(itertools.chain(*map(split_segment, segments)))

def clip_to_window(postings, lon_center):
    tx = lambda p: (p[0], geodesy.anglenorm(p[1], 180. - lon_center))
    segment = [tx(postings[0])]
    for a, b in pairwise_postings(postings):
        last_p = segment[-1]
        p = tx(b)
        if abs(p[1] - last_p[1]) > 180.:
            dir = (-1 if last_p[1] < lon_center else 1)
            cross_in = lon_center + 180 * dir
            cross_out = lon_center - 180 * dir
            midlat = last_p[0] + (p[0] - last_p[0]) * (cross_in - last_p[1]) / (geodesy.anglenorm(p[1], 180 - last_p[1]) - last_p[1])
            segment.append((midlat, cross_in))
            yield segment
            segment = [(midlat, cross_out)]
        segment.append(p)
    yield segment

def mercator_safe_linear_segment(start, end, tolerance):
    def midpoint(a, b):
        xyza = geodesy.ll_to_ecefu(a)
        xyzb = geodesy.ll_to_ecefu(b)
        xyzmid = geodesy.vnorm(geodesy.vadd(xyza, xyzb))
        return geodesy.ecefu_to_ll(xyzmid)

    def merc_midpoint(a, b):
        merca = math.log(math.tan(math.pi/4 + math.radians(a[0])/2.))
        mercb = math.log(math.tan(math.pi/4 + math.radians(b[0])/2.))
        mercmid = .5*(merca + mercb)
        midlat = math.degrees(2*math.atan(math.exp(mercmid)) - math.pi/2)
        midlon = geodesy.anglenorm(a[1] + .5*geodesy.anglenorm(b[1] - a[1]))
        return (midlat, midlon)
    
    MIN_SEG_LENGTH = 1. # m

    def interim_points(start, end):
        if geodesy.distance(start, end) < MIN_SEG_LENGTH:
            return

        mid = midpoint(start, end)
        merc_mid = merc_midpoint(start, end)
        if geodesy.distance(mid, merc_mid) <= tolerance:
            return

        for p in interim_points(start, mid):
            yield p
        yield mid
        for p in interim_points(mid, end):
            yield p
        
    yield start
    for p in interim_points(start, end):
        yield p
    yield end

def make_mercator_safe(postings, tolerance):
    yield postings[0]
    for a, b in pairwise_postings(postings):
        for p in list(mercator_safe_linear_segment(a, b, tolerance))[1:]:
            yield p

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
        (r'/geojson/(?P<tag>.*)', GeojsonHandler, {}, 'geojson'),
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
