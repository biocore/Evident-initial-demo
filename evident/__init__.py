#!/usr/bin/env python

__author__ = "e-vident Development Team"
__copyright__ = "Copyright 2011-2012, e-vident Project"
__credits__ = ["Antonio Gonzalez Pena", "Meg Pirrung", "Will Van Treuren"]
__license__ = "GPL"
__version__ = "0.9"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

__all__ = ['subsampling','map_sample_space']

import logging
from os import chmod, environ

# if level is set to DEBUG log messages will be written
logging.basicConfig(filename='/tmp/e-vident.log', level=logging.DEBUG, \
                    format='[%(asctime)s].%(levelname)s: %(message)s')

# when the file is created assure everyone can w/r to it 
try:
    chmod('/tmp/e-vident.log',0777)
except:
    pass

# performance get's affected by the biom backend used, forced to use SparseMat
environ['BIOM_CONFIG_FP'] = '/evident/data/biom_config'