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
import alt
import time

from munsell import munsell as m
m.init()

IX = None

class LandfallHandler(web.RequestHandler):
    def get(self): # should be POST
        _origin = self.get_argument('origin')
        size = int(self.get_argument('size'))
        _range = self.get_argument('range', '0,360')
        mindist = float(self.get_argument('mindist', '10000'))

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

        global IX
        if IX is None:
            print 'loading index...'
            IX = pickle.load(open('data/tmp/tagged_coastline'))

        print 'generating...'
        postings = alt.postings(origin, IX, res, range[0], range[1], mindist)

        print 'saving...'
        output = {
            'origin': origin,
            'res': res,
            'range': range,
            'min_dist': mindist,
            'postings': [[dist, list(areas)] for dist, areas in postings],
        }
        tag = '%d.json' % time.time()
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
        self.render('render.html', data=data)

class KmlHandler(web.RequestHandler):
    def get(self, tag):
        with open(os.path.join(config.OUTPUT_PATH, tag)) as f:
            data = json.load(f)

        antipode = (-data['origin'][0], geodesy.anglenorm(data['origin'][1] + 180.))

        path = []
        prevdist = None
        for i, (dist, _) in enumerate(itertools.chain(data['postings'], [data['postings'][0]])):
            if dist < 0:
                dist = alt.EARTH_CIRCUMF - data['min_dist'] # TODO go all the way back to start but insert interstitial point
            bear = data['range'][0] + i * data['res']
            p = geodesy.plot(data['origin'], bear, dist)[0]
            if prevdist is not None and (dist > .5*alt.EARTH_CIRCUMF) != (prevdist > .5*alt.EARTH_CIRCUMF):
                path.append(antipode)
            path.append(p)
            prevdist = dist

        SEGSIZE = 1000
        segments = [path[i:min(i+1+SEGSIZE, len(path))] for i in xrange(0, len(path), SEGSIZE)]
        segments.append([geodesy.plot(data['origin'], bear, data['min_dist'])[0] for bear in xrange(0, 361, 5)])

        self.set_header('Content-Type', 'application/vnd.google-earth.kml+xml')
        self.set_header('Content-Disposition', 'attachment; filename="landfall.kml"')
        self.render('render.kml', segments=segments)

class ColorProviderHandler(websocket.WebSocketHandler):
    def open(self):
        pass

    def check_origin(self, origin):
        return True

    def on_message(self, message):
        data = json.loads(message)

        num_colors = data['num_colors']
        LUM = data['lum']
        SATFAR, SATCLOSE = data['sat']
        STEPS = data['steps']

        def color(h, v, c):
            rgb = m.munsell(h, v, c)
            if not m.in_gamut(rgb):
                return None
            return '#%s' % ''.join('%02x' % k for k in m.rgb_to_hex(rgb))

        hues = [float(i)/num_colors for i in range(num_colors)]
        steps = [float(i)/(STEPS-1) for i in range(STEPS)]
        table = [[color(hue, LUM, SATCLOSE*(1-k)+SATFAR*k) for k in steps] for hue in hues]

        self.write_message(json.dumps(table))

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
        (r'/list', OutputListHandler),
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
