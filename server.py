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
        tag = '%s.json' % datetime.now().strftime('%Y%m%d%H%M%S')
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
        with open(os.path.join(config.OUTPUT_PATH, tag)) as f:
            data = json.load(f)
        self.render('render.html', data=data, info=self.get_admin_info(data))

    def get_admin_info(self, data):
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

        entities = set((posting[1] or pd.CC_UNCLAIMED) for posting in data['postings'])
        return dict((k, v) for k, v in info.iteritems() if k in entities)

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

class ColorProviderHandler(websocket.WebSocketHandler):
    def open(self):
        pass

    def check_origin(self, origin):
        return True

    def on_message(self, message):
        data = json.loads(message)

        hues = data['hues']
        lummin, lummax = data['lum']
        satmin, satmax = data['sat']

        def get_ramp(hue):
            key = (hue, lummin, lummax, satmin, satmax)
            if key not in COLOR_CACHE:
                print 'generating color ramp: hue %s lum %s %s sat %s %s' % key
                COLOR_CACHE[key] = color_ramp(*key)
            return COLOR_CACHE[key]

        colors = dict((hue, get_ramp(hue)) for hue in hues)
        self.write_message(json.dumps(colors))

    def on_close(self):
        pass

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
        (r'/colors', ColorProviderHandler),
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
