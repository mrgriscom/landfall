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

from munsell import munsell as m
m.init()

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
