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
from evident.rarefaction import _format_rarefactions

class TopLevelTests(TestCase):
    def setUp(self):
        pass

    def test_format_rarefactions(self):
        sample_ids = ['M9Thmr217.141049', 'M3Midl217.141103', \
            'M9Midr217.141060', 'M3Pinr217.141057']
        output = _format_rarefactions(metrics_data_a, sample_ids)
        self.assertEquals(output, expected_formated_dict_a)
        pass

metrics_data_a = {\
'PD_whole_tree': [\
    ['alpha_rarefaction_381_0', 381, 0.0, '6.62714', '6.71907', '6.03955', '7.18497'],\
    ['alpha_rarefaction_666_0', 666, 0.0, '7.41714', '7.44263', '6.26622', '9.85806'],\
    ['alpha_rarefaction_951_0', 951, 0.0, '9.72276', '8.49664', '7.27135', '10.89548'],\
    ['alpha_rarefaction_1236_0', 1236, 0.0, '9.43287', '9.88928', '7.68549', '11.97569'],\
    ['alpha_rarefaction_1521_0', 1521, 0.0, '10.14783', '10.28997', '8.04553', '12.79367']\
],\
'observed_species': [\
    ['alpha_rarefaction_381_0', 381, 0, '49.0', '45.0', '39.0', '66.0'],\
    ['alpha_rarefaction_666_0', 666, 0, '63.0', '66.0', '42.0', '94.0'],\
    ['alpha_rarefaction_951_0', 951, 0, '79.0', '74.0', '56.0', '118.0'],\
    ['alpha_rarefaction_1236_0', 1236, 0, '83.0', '84.0', '60.0', '140.0'],\
    ['alpha_rarefaction_1521_0', 1521, 0, '90.0', '100.0', '70.0', '150.0']\
],\
'chao1': [\
    ['alpha_rarefaction_381_0', 381, 0, '99.1428571429', '90.1111111111', '114.6', '160.6'],\
    ['alpha_rarefaction_666_0', 666, 0, '117.090909091', '123.0', '88.4285714286', '153.913043478'],\
    ['alpha_rarefaction_951_0', 951, 0, '146.363636364', '184.0', '81.2', '233.0'],\
    ['alpha_rarefaction_1236_0', 1236, 0, '104.565217391', '174.461538462', '80.6470588235', '261.538461538'],\
    ['alpha_rarefaction_1521_0', 1521, 0, '121.166666667', '192.8125', '95.8333333333', '283.956521739']\
]}

expected_formated_dict_a={\
'PD_whole_tree':\
    (['', 'sequences per sample', 'iteration', 'M9Thmr217.141049', 'M3Midl217.141103', 'M9Midr217.141060', 'M3Pinr217.141057'], [],\
    ['alpha_rarefaction_381_0', 'alpha_rarefaction_666_0', 'alpha_rarefaction_951_0', 'alpha_rarefaction_1236_0', 'alpha_rarefaction_1521_0'],\
    [\
    [381.0, 0.0, 6.62714, 6.71907, 6.03955, 7.18497],\
    [666.0, 0.0, 7.41714, 7.44263, 6.26622, 9.85806],\
    [951.0, 0.0, 9.72276, 8.49664, 7.27135, 10.89548],\
    [1236.0, 0.0, 9.43287, 9.88928, 7.68549, 11.97569],\
    [1521.0, 0.0, 10.14783, 10.28997, 8.04553, 12.79367]]),\
'observed_species':\
    (['', 'sequences per sample', 'iteration', 'M9Thmr217.141049', 'M3Midl217.141103', 'M9Midr217.141060', 'M3Pinr217.141057'], [],\
    ['alpha_rarefaction_381_0', 'alpha_rarefaction_666_0', 'alpha_rarefaction_951_0', 'alpha_rarefaction_1236_0', 'alpha_rarefaction_1521_0'],\
    [[381.0, 0.0, 49.0, 45.0, 39.0, 66.0],\
    [666.0, 0.0, 63.0, 66.0, 42.0, 94.0],\
    [951.0, 0.0, 79.0, 74.0, 56.0, 118.0],\
    [1236.0, 0.0, 83.0, 84.0, 60.0, 140.0],\
    [1521.0, 0.0, 90.0, 100.0, 70.0, 150.0]]),\
'chao1':\
    (['', 'sequences per sample', 'iteration', 'M9Thmr217.141049', 'M3Midl217.141103', 'M9Midr217.141060', 'M3Pinr217.141057'], [],\
    ['alpha_rarefaction_381_0', 'alpha_rarefaction_666_0', 'alpha_rarefaction_951_0', 'alpha_rarefaction_1236_0', 'alpha_rarefaction_1521_0'],\
    [[381.0, 0.0, 99.1428571429, 90.1111111111, 114.6, 160.6],\
    [666.0, 0.0, 117.090909091, 123.0, 88.4285714286, 153.913043478],\
    [951.0, 0.0, 146.363636364, 184.0, 81.2, 233.0],\
    [1236.0, 0.0, 104.565217391, 174.461538462, 80.6470588235, 261.538461538],\
    [1521.0, 0.0, 121.166666667, 192.8125, 95.8333333333, 283.956521739]])\
}


if __name__ == "__main__":
    main()
