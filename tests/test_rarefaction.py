#!/usr/bin/env python
# File created on 26 Sep 2012
from __future__ import division

__author__ = "Yoshiki Vazquez-Baeza"
__copyright__ = "Copyright 2011-2012, Evident"
__credits__ = ["Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = "0.9-dev"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "yoshiki89@gmail.com"
__status__ = "Development"


from cogent.util.unit_test import TestCase, main

from shutil import rmtree
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, create_dir
from qiime.util import (get_qiime_temp_dir, 
                        get_tmp_filename)
from qiime.test import initiate_timeout, disable_timeout

class TopLevelTests(TestCase):
    
    def setUp(self):
        
        self.files_to_remove = []
        self.dirs_to_remove = []
        
        # Create example output directory
        tmp_dir = get_qiime_temp_dir()
        self.test_out = get_tmp_filename(tmp_dir=tmp_dir,
                                         prefix='qiime_parallel_tests_',
                                         suffix='',
                                         result_constructor=str)
        self.dirs_to_remove.append(self.test_out)
        create_dir(self.test_out)
        
        # Create example input file
        self.inseqs1_fp = get_tmp_filename(tmp_dir=self.test_out,
                                            prefix='qiime_inseqs',
                                            suffix='.fasta')
        inseqs1_f = open(self.inseqs1_fp,'w')
        inseqs1_f.write(inseqs1)
        inseqs1_f.close()
        self.files_to_remove.append(self.inseqs1_fp)
        
        # Define number of seconds a test can run for before timing out 
        # and failing
        initiate_timeout(60)

    
    def tearDown(self):
        
        disable_timeout()
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

inseqs1 = """>example input here
ACGT
"""




if __name__ == "__main__":
    main()
