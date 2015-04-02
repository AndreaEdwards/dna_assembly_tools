#!/usr/bin/env python

__author__ = "Andrea L Halweg-Edwards"
__copyright__ = "Copyright 2015, The LASER Project"
__credits__ = ["Andrea L Halweg-Edwards"]
__license__ = "BSD"
__version__ = "0.1.0-dev"
__maintainer__ = "Andrea L Halweg-Edwards"
__email__ = "andrea.edwards@colorado.edu"
__status__ = "Development"


class ConfigurationManager(object):
    """"""
    def __init__(self):
        self.web_port = 3005

config = ConfigurationManager()
