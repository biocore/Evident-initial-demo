#!/usr/bin/env python
# File created on 12 Nov 2012
#from __future__ import division

__author__ = "Antonio Gonzalez"
__copyright__ = "Copyright 2012, Evident"
__credits__ = ["Antonio Gonzalez"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "Antonio Gonzalez"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

"""code for doing error handling"""

from mod_python import Session

def raiseApacheError(text):
	"""Receives a text field and raises an exception in apache"""
	req.write(text)
	raise apache.SERVER_RETURN, 0
