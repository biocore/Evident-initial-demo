#!/usr/bin/env python
# File created on 29 Jan 2013
from __future__ import division

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2013, The Evident Project"
__credits__ = ["Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"


from shutil import rmtree
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, create_dir
from qiime.util import (get_qiime_temp_dir, 
                        get_tmp_filename)
from qiime.test import initiate_timeout, disable_timeout
from numpy import isnan, nan, array, std, allclose
from evident.compare_treatment_dists import (between_treatments_dist, 
    treatment_dist, treatment_covering, within_treatment_dist, 
    compare_treatment_dists)
from cogent.parse.tree import DndParser
from biom.parse import parse_biom_table_str

class compare_treatments(TestCase):
    
    def setUp(self):
        """Create data used by the majority of tests."""
        self.distmat = array([[ 0.        ,  0.00071582,  0.00071582,  0.05726557,  0.05010737,
         0.07158196,  0.16535433,  0.07874016,  0.07874016,  0.08231926,
         0.10021475,  0.02219041,  0.00071582],
       [ 0.00071582,  0.        ,  0.        ,  0.05798139,  0.05082319,
         0.07229778,  0.16607015,  0.07945598,  0.07945598,  0.08303508,
         0.10093057,  0.02148997,  0.        ],
       [ 0.00071582,  0.        ,  0.        ,  0.05798139,  0.05082319,
         0.07229778,  0.16607015,  0.07945598,  0.07945598,  0.08303508,
         0.10093057,  0.02148997,  0.        ],
       [ 0.05726557,  0.05798139,  0.05798139,  0.        ,  0.02243829,
         0.04487659,  0.14285714,  0.13600573,  0.13600573,  0.13958482,
         0.15748031,  0.07945598,  0.05798139],
       [ 0.05010737,  0.05082319,  0.05082319,  0.02243829,  0.        ,
         0.02260739,  0.1213263 ,  0.11535689,  0.11535689,  0.11896179,
         0.1369863 ,  0.07229778,  0.05082319],
       [ 0.07158196,  0.07229778,  0.07229778,  0.04487659,  0.02260739,
         0.        ,  0.10100231,  0.1369863 ,  0.1369863 ,  0.1405912 ,
         0.11790715,  0.05193855,  0.07229778],
       [ 0.16535433,  0.16607015,  0.16607015,  0.14285714,  0.1213263 ,
         0.10100231,  0.        ,  0.09401709,  0.09401709,  0.0979021 ,
         0.07239459,  0.14776884,  0.16607015],
       [ 0.07874016,  0.07945598,  0.07945598,  0.13600573,  0.11535689,
         0.1369863 ,  0.09401709,  0.        ,  0.        ,  0.003885  ,
         0.02331002,  0.10093057,  0.07945598],
       [ 0.07874016,  0.07945598,  0.07945598,  0.13600573,  0.11535689,
         0.1369863 ,  0.09401709,  0.        ,  0.        ,  0.003885  ,
         0.02331002,  0.10093057,  0.07945598],
       [ 0.08231926,  0.08303508,  0.08303508,  0.13958482,  0.11896179,
         0.1405912 ,  0.0979021 ,  0.003885  ,  0.003885  ,  0.        ,
         0.02719503,  0.10450966,  0.08303508],
       [ 0.10021475,  0.10093057,  0.10093057,  0.15748031,  0.1369863 ,
         0.11790715,  0.07239459,  0.02331002,  0.02331002,  0.02719503,
         0.        ,  0.08119971,  0.10093057],
       [ 0.02219041,  0.02148997,  0.02148997,  0.07945598,  0.07229778,
         0.05193855,  0.14776884,  0.10093057,  0.10093057,  0.10450966,
         0.08119971,  0.        ,  0.02148997],
       [ 0.00071582,  0.        ,  0.        ,  0.05798139,  0.05082319,
         0.07229778,  0.16607015,  0.07945598,  0.07945598,  0.08303508,
         0.10093057,  0.02148997,  0.        ]])
        self.samples = ['a1', 'a2', 'a3', 'b1', 'b2', 'b3', 'c1', 'c2', 'c3', 'c4', 
            'd1', 'd2', 'd3']

    def test_between_treatments_dist(self):
        """Tests that distance between treatment groups is calculated right."""
        g1 = ['a3','a1','b3']
        g2 = ['c1', 'd2', 'c3']
        exp_m = 0.091469795660200012
        exp_se = 0.0058270292521978868
        obs_m, obs_se = between_treatments_dist(g1, g2, self.samples, 
            self.distmat)
        self.assertFloatEqual(exp_m, obs_m)
        self.assertFloatEqual(exp_se, obs_se)

    def test_treatment_dist(self):
        """Tests that intersample distance is calculated correctly."""
        group = ['a1','a3','b3']
        exp_m = 0.072436576333380465
        exp_se = 0.0015658874661095082
        obs_m, obs_se = treatment_dist(group, self.samples, self.distmat)
        self.assertFloatEqual(exp_m, obs_m)
        self.assertFloatEqual(exp_se, obs_se)
        # test if ordering is different
        group = ['a3','a1','b3']
        obs_m, obs_se = treatment_dist(group, self.samples, self.distmat)
        self.assertFloatEqual(exp_m, obs_m)
        self.assertFloatEqual(exp_se, obs_se)
        # test when group=samples, shoudl return nan nan
        group = self.samples
        obs_m, obs_se = treatment_dist(group, self.samples, self.distmat)
        self.assertTrue(isnan(obs_m))
        self.assertTrue(isnan(obs_se))
        # test with 5 samples
        group = ['a3', 'b1', 'b2','d3','a1']
        exp_m = 0.086332301009816426
        exp_se = 0.0011588621239506576
        obs_m, obs_se = treatment_dist(group, self.samples, self.distmat)
        self.assertFloatEqual(exp_m, obs_m)
        self.assertFloatEqual(exp_se, obs_se)

    def test_within_treatment_dist(self):
        """Tests that intrasample distance is calculated correctly."""
        # test one group in different order
        group1 = ['a1','a2','c1']
        obs_m, obs_se = within_treatment_dist(group1, self.samples, self.distmat)
        exp_m = 0.11071343333333333
        exp_se = 0.012963434547075115
        self.assertFloatEqual(obs_m, exp_m)
        self.assertFloatEqual(obs_se, exp_se)
        group1 = ['a1','c1','a2'] #show order doesnt matter
        obs_m, obs_se = within_treatment_dist(group1, self.samples, self.distmat)
        self.assertFloatEqual(obs_m, exp_m)
        self.assertFloatEqual(obs_se, exp_se)
        # test that nans returned when the se is masked and the mean is nan
        # this happens when you get treatments represented by only a single 
        # sample
        group1 = ['d1']
        obs_m, obs_se = within_treatment_dist(group1, self.samples, self.distmat)
        self.assertTrue(isnan(obs_m))
        self.assertTrue(isnan(obs_se))

    def test_treatment_covering(self):
        """Tests treatment covering returns the correct data."""
        #sids = [['a1', 'a2'], ['c1'], ['d1', 'd2', 'd3']]
        sids = ['a1', 'a2', 'c1', 'd1', 'd2', 'd3']
        
        mf = \
            {'a1': {'Diet': 'LF', 'HSID': 'a', 'Pref': '1'},
             'a2': {'Diet': 'LF', 'HSID': 'a', 'Pref': '1'},
             'a3': {'Diet': 'HF', 'HSID': 'a', 'Pref': '1'},
             'b1': {'Diet': 'HF', 'HSID': 'b', 'Pref': '4'},
             'b2': {'Diet': 'HF', 'HSID': 'b', 'Pref': '5'},
             'b3': {'Diet': 'HF', 'HSID': 'b', 'Pref': '5'},
             'c1': {'Diet': 'LF', 'HSID': 'c', 'Pref': '5'},
             'c2': {'Diet': 'LF', 'HSID': 'c', 'Pref': '5'},
             'c3': {'Diet': 'LF', 'HSID': 'c', 'Pref': '5'},
             'c4': {'Diet': 'LF', 'HSID': 'd', 'Pref': '2'},
             'd1': {'Diet': 'HF', 'HSID': 'd', 'Pref': '2'},
             'd2': {'Diet': 'HF', 'HSID': 'd', 'Pref': '1'},
             'd3': {'Diet': 'HF', 'HSID': 'd', 'Pref': '1'}}
        exp = {'HF': ['d1', 'd2', 'd3'], 'LF': ['a1', 'a2', 'c1']}
        self.assertEqual(exp, treatment_covering(sids, 'Diet', mf))
        exp = {'a': ['a1', 'a2'], 'c': ['c1'], 'd': ['d1', 'd2', 'd3']}
        self.assertEqual(exp, treatment_covering(sids, 'HSID', mf))
        exp = {'1': ['a1', 'a2', 'd2', 'd3'], '2': ['d1'], '5': ['c1']}
        self.assertEqual(exp , treatment_covering(sids, 'Pref', mf))
        # test with a single sample
        sids = ['a1']
        exp = {'1': ['a1']}
        self.assertEqual(exp , treatment_covering(sids, 'Pref', mf))

    def test_compare_treatment_dists(self):
        """Tests the the entire library functions as expected."""
        bt = parse_biom_table_str('{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "test_code","date": "2013-01-30T14:47:00.868638","matrix_type": "sparse","matrix_element_type": "float","shape": [8, 13],"data": [[0,0,100.0],[0,1,12.0],[0,2,45.0],[0,7,24.0],[0,8,67.0],[0,9,132.0],[0,10,991.0],[0,11,21.0],[0,12,5.0],[1,0,54.0],[1,1,989.0],[1,2,200.0],[1,3,425.0],[1,4,2.0],[1,5,4.0],[1,11,52.0],[1,12,13.0],[2,0,4.0],[2,3,11.0],[2,4,11.0],[2,5,100.0],[2,6,491.0],[2,7,55.0],[2,8,98.0],[2,9,54.0],[2,10,104.0],[3,0,45.0],[3,1,78.0],[3,2,52.0],[3,3,14.0],[3,11,1000.0],[3,12,141.0],[4,0,92.0],[4,1,46.0],[4,2,49.0],[4,3,770.0],[4,4,14.0],[4,7,1.0],[4,8,55.0],[4,9,255.0],[4,12,14.0],[5,0,67.0],[5,1,92.0],[5,2,800.0],[5,4,13.0],[5,5,17.0],[5,6,29.0],[5,7,27.0],[5,8,25.0],[5,9,221.0],[5,10,228.0],[5,11,9.0],[5,12,9.0],[6,0,150.0],[6,1,149.0],[6,2,11.0],[6,3,35.0],[6,4,899.0],[6,5,766.0],[6,6,348.0],[6,7,680.0],[6,8,496.0],[6,9,467.0],[6,10,13.0],[6,11,327.0],[6,12,855.0],[7,0,300.0],[7,1,45.0],[7,2,78.0],[7,3,22.0],[7,4,361.0],[7,5,271.0],[7,6,531.0],[7,7,256.0],[7,8,251.0],[7,10,70.0],[7,11,250.0],[7,12,31.0]],"rows": [{"id": "O1", "metadata": null},{"id": "O2", "metadata": null},{"id": "O3", "metadata": null},{"id": "O4", "metadata": null},{"id": "O5", "metadata": null},{"id": "O6", "metadata": null},{"id": "O7", "metadata": null},{"id": "O8", "metadata": null}],"columns": [{"id": "a1", "metadata": null},{"id": "a2", "metadata": null},{"id": "a3", "metadata": null},{"id": "b1", "metadata": null},{"id": "b2", "metadata": null},{"id": "b3", "metadata": null},{"id": "c1", "metadata": null},{"id": "c2", "metadata": null},{"id": "c3", "metadata": null},{"id": "c4", "metadata": null},{"id": "d1", "metadata": null},{"id": "d2", "metadata": null},{"id": "d3", "metadata": null}]}')
        mf = \
            {'a1': {'Diet': 'LF', 'HSID': 'a', 'Pref': '1'},
             'a2': {'Diet': 'LF', 'HSID': 'a', 'Pref': '1'},
             'a3': {'Diet': 'HF', 'HSID': 'a', 'Pref': '1'},
             'b1': {'Diet': 'HF', 'HSID': 'b', 'Pref': '4'},
             'b2': {'Diet': 'HF', 'HSID': 'b', 'Pref': '5'},
             'b3': {'Diet': 'HF', 'HSID': 'b', 'Pref': '5'},
             'c1': {'Diet': 'LF', 'HSID': 'c', 'Pref': '5'},
             'c2': {'Diet': 'LF', 'HSID': 'c', 'Pref': '5'},
             'c3': {'Diet': 'LF', 'HSID': 'c', 'Pref': '5'},
             'c4': {'Diet': 'LF', 'HSID': 'd', 'Pref': '2'},
             'd1': {'Diet': 'HF', 'HSID': 'd', 'Pref': '2'},
             'd2': {'Diet': 'HF', 'HSID': 'd', 'Pref': '1'},
             'd3': {'Diet': 'HF', 'HSID': 'd', 'Pref': '1'}}
        tr = DndParser('(((O1:0.06,O2:0.1)A:0.031,(O3:0.001,O4:0.01)B:0.2)AB:0.4,((O5:0.03,O6:0.02)C:0.13,(O7:0.01,O8:0.005)D:0.1)CD:0.3)root;')
        # test known quantitiy 
        sids = [['a1', 'a2'], ['c1'], ['d1', 'd2', 'd3']]
        sids = ['a1', 'a2','c1', 'd1', 'd2', 'd3']
        
        exp_out = (['LF', 'HF'],
            array([[ 0.11071343,  0.07019723],
               [ 0.        ,  0.06787341]]),
            array([[ 0.01296343,  0.00657256],
               [ 0.        ,  0.00562879]]),
            array([[ 0.0763829 ,  0.0015125 ],
               [ 0.0724499 ,  0.00154609]]))
        obs_out = compare_treatment_dists(sids, 'Diet', mf, bt, 
            'unweighted_unifrac', tr)
        self.assertEqual(exp_out[0], obs_out[0])
        for i,j in zip(obs_out[1:], exp_out[1:]):
            self.assertTrue(allclose(i,j)) # have to use all close bc tiny float errors

if __name__ == "__main__":
    main()