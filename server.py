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
import bisect

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

class RenderHandler(web.RequestHandler):
    def get(self, tag):
        width = int(self.get_argument('width', '3000'))
        height = int(self.get_argument('height', '800'))

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

        force_interfere = self.get_argument('diffcolor', '')
        force_interfere = [pair.split(':') for pair in force_interfere.split(',') if pair]
        no_subdivisions = self.get_argument('nosubdiv', '')
        no_subdivisions = set(filter(None, no_subdivisions.split(',')))
        resolve_as = self.get_argument('resolve', '')
        resolve_as = dict(pair.split(':') for pair in resolve_as.split(',') if pair)

        dist_unit = self.get_argument('distunit', 'km')
        assert dist_unit in ('km', 'mi', 'deg', 'none')

        params = {
            'dim': [width, height],
            'colors': colors,
            'force_interfere': force_interfere,
            'no_subdivisions': no_subdivisions,
            'resolve_as': resolve_as,
            'dist_unit': dist_unit,
        }

        with open(os.path.join(config.OUTPUT_PATH, tag)) as f:
            data = json.load(f)
        data['size'] = len(data['postings'])
        data['lonspan'] = (data['range'][1] - data['range'][0]) % 360. or 360.
        data['wraparound'] = (data['lonspan'] == 360.)
        self.process_admins(data, params)

        self.render('render.html', data=data, params=params)

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

        def change_admin(a):
            a = params['resolve_as'].get(a, a)
            parent = admin_info[a]['parent']
            if parent and parent in params['no_subdivisions']:
                a = parent
            return a
        admin_postings = map(change_admin, admin_postings)
        admins = set(admin_postings)
        pprint(dict((a, admin_info[a]) for a in admins))

        def ix_offset(i, offset, size, wrap=data['wraparound']):
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
            assert a in admins and b in admins
            interfere(a, b)

        MAX_ADJ_INTERF = 2
        for delta in xrange(1, MAX_ADJ_INTERF + 1):
            interf_prio.append('adj%d' % delta)
            for i in xrange(len(admin_segments)):
                admin = admin_segments[i]['admin']
                for dir in (-1, 1):
                    adj = ix_offset(i, dir*delta, len(admin_segments))
                    if adj < 0:
                        continue
                    adj_admin = admin_segments[adj]['admin']
                    interfere(admin, adj_admin)

        MAX_PX_INTERF = 50
        interf_prio.append('nearpx')
        unwrapped = list(map(dict, admin_segments)) # deep copy
        if unwrapped:
            first = unwrapped[0]
            if first['start'] > first['end']:
                unwrapped.append({'admin': first['admin'], 'start': first['start'], 'end': data['size'] - 1})
                first['start'] = 0
        seg_starts = [seg['start'] for seg in unwrapped]
        seg_ends = [seg['end'] for seg in unwrapped]
        def overlapping_segments(min, max):
            assert 0 <= min <= max < data['size']
            start = bisect.bisect_left(seg_ends, min)
            end = bisect.bisect_right(seg_starts, max) - 1
            return unwrapped[start:end+1]
        world_width = int(round(360. / data['res']))
        assert world_width >= data['size']
        assert not data['wraparound'] or world_width == data['size']
        for i, seg in enumerate(admin_segments):
            seg_width = (seg['end'] - seg['start']) % world_width + 1
            if seg_width + 2*MAX_PX_INTERF >= world_width:
                ranges = [(0, world_width - 1)]
            else:
                min = (seg['start'] - MAX_PX_INTERF) % world_width
                max = (seg['end'] + MAX_PX_INTERF) % world_width
                if min <= max:
                    ranges = [(min, max)]
                else:
                    ranges = [(0, max), (min, world_width - 1)]
            for ovl in itertools.chain(*(overlapping_segments(*rng) for rng in ranges)):
                interfere(seg['admin'], ovl['admin'])

        pprint(interferences)




        """


    // build adjacency graph
    var adj = {};
    var adjedge = function(a, b) {
        if (!adj[a]) {
            adj[a] = {};
        }
        adj[a][b] = true;
    };
    for (var i = 0; i < DATA.size; i++) {
        var a0 = admin_postings[i];
        var a1 = admin_postings[(i + 1) % admin_postings.length];
        if (a0 == a1) {
            continue;
        }
        adjedge(a0, a1);
        adjedge(a1, a0);
    }
    _.each(adj, function(v, k) {
        adj[k] = _.keys(v);
    });

    var colors = assign_colors(PARAMS.hues.length, admins, adj);
    console.log(colors);
}
"""




class KmlHandler(web.RequestHandler):
    def get(self, tag):
        with open(os.path.join(config.OUTPUT_PATH, tag)) as f:
            data = json.load(f)

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
