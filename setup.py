#!/usr/bin/env python

from distutils.core import setup
from glob import glob

__author__ = "Yoshiki Vazquez Baeza"
__copyright__ = "Copyright 2013, The Evident Project"
__credits__ = ["Antonio Gonzalez Pena", "Meg Pirrung", "William Van Treuren",
               "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "Yoshiki Vazquez Baeza"
__email__ = "yoshiki89@gmail.com"
__status__ = "Development"

# based on the text found in github.com/qiime/pynast
classes = """
    Development Status :: 4
    Beta License :: OSI Approved :: GPL License
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: User Interfaces
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: OS Independent
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""

classifiers = [s.strip() for s in classes.split('\n') if s]

long_description = ("Evident: elucidating sampling effort for microbial "
                    "analysis studies")

setup(name='evident',
        version=__version__,
        description='Evident',
        author=("Antonio Gonzalez Pena, Meg Pirrung, William Van Treuren & "
                "Yoshiki Vazquez Baeza"),
        author_email=__email__,
        maintainer=__maintainer__,
        maintainer_email=__email__,
        url='http://github.com/biocore/evident',
        packages=['evident'],
        scripts=glob('scripts/*py'),
        package_data={},
        data_files=('data',['biom_config',
                    'bushman_enterotypes_cafe_alpha_stddev.html',
                    'bushman_enterotypes_cafe_alpha_stderr.html',
                    'bushman_enterotypes_cafe.biom',
                    'bushman_enterotypes_cafe_map.txt',
                    'bushman_enterotypes_cafe_selectors.txt',
                    'bushman_enterotypes_cafe_unweighted_unifrac_pc.txt',
                    'bushman_enterotypes_cafe_weighted_unifrac_pc.txt',
                    'gg_97_otus_4feb2011.tre',
                    'hmp-v13_alpha_stddev.html',
                    'hmp-v13_alpha_stderr.html',
                    'hmp-v13.biom',
                    'hmp-v13_map.txt',
                    'hmp-v13_selectors.txt',
                    'hmp-v13_unweighted_unifrac_pc.txt',
                    'hmp-v35_alpha_stddev.html',
                    'hmp-v35_alpha_stderr.html',
                    'hmp-v35.biom',
                    'hmp-v35_map.txt',
                    'hmp-v35_selectors.txt',
                    'hmp-v35_unweighted_unifrac_pc.txt',
                    'keyboard_alpha_stddev.html',
                    'keyboard_alpha_stderr.html',
                    'keyboard.biom',
                    'keyboard_map.txt',
                    'keyboard_selectors.txt',
                    'keyboard_unweighted_unifrac_pc.txt',
                    'studies.txt',
                    'turnbaugh-twins_alpha_stddev.html',
                    'turnbaugh-twins_alpha_stderr.html',
                    'turnbaugh-twins.biom',
                    'turnbaugh-twins_map.txt',
                    'turnbaugh-twins_selectors.txt',
                    'turnbaugh-twins_unweighted_unifrac_pc.txt',
                    'wholebody_alpha_stddev.html',
                    'wholebody_alpha_stderr.html',
                    'wholebody.biom',
                    'wholebody_map.txt',
                    'wholebody_selectors.txt',
                    'wholebody_unweighted_unifrac_pc.txt']),
        install_requires=["numpy >= 1.5.1, <=1.7.1", 'qiime'],
        long_description=long_description,
        classifiers=classifiers)

