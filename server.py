#!/usr/bin/env python

from os.path import dirname, join
from base64 import b64encode
from uuid import uuid4

import tornado.web
from tornado.ioloop import IOLoop
from tornado.httpserver import HTTPServer
from tornado.options import define, options, parse_command_line

from config import *
from handlers.base_handlers import (
    IndexHandler
    )
from handlers.tool_handlers import (
    CPECHandler, UploadHandler
    )


class Application(tornado.web.Application):
    def __init__(self):
        DIRNAME = dirname(__file__)
        STATIC_PATH = join(DIRNAME, "static")
        TEMPLATE_PATH = join(DIRNAME, "templates")  # base folder for webpages
        RES_PATH = join(DIRNAME, "results")
        COOKIE_SECRET = b64encode(uuid4().bytes + uuid4().bytes)
        handlers = [
                (r"/", IndexHandler),
                (r"/results/(.*)", tornado.web.StaticFileHandler,
                 {"path": RES_PATH}),
                (r"/static/(.*)", tornado.web.StaticFileHandler,
                 {"path": STATIC_PATH}),
                (r"/assembly_designer/", CPECHandler),
                (r"/upload/", UploadHandler),
                ]
        settings = {
            "template_path": TEMPLATE_PATH,
            "debug": True,
            "cookie_secret": COOKIE_SECRET,
            "login_url": "/",
        }
        tornado.web.Application.__init__(self, handlers, **settings)


def main():
    define("port", default=config.web_port, type=int)
    parse_command_line()
    http_server = HTTPServer(Application())
    http_server.listen(options.port)
    print "Tornado started on port", options.port
    IOLoop.instance().start()

if __name__ == "__main__":
    main()
