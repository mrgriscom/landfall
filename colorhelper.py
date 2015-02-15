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

import munsell as m
m.init()

class WebSocketTestHandler(websocket.WebSocketHandler):
    def open(self):
        pass

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

    try:
        port = int(sys.argv[1])
    except IndexError:
        port = 8000

    application = web.Application([
        (r'/socket', WebSocketTestHandler),
    ])
    application.listen(port)

    try:
        IOLoop.instance().start()
    except KeyboardInterrupt:
        pass
    except Exception, e:
        print e
        raise

    logging.info('shutting down...')
