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
from qiime.util import get_qiime_temp_dir, get_tmp_filename

from biom.parse import parse_biom_table
from biom.exception import TableException
from evident.subsampling import select_samples

class TopLevelTests(TestCase):
    
    def setUp(self):
        self.biom_object = parse_biom_table(input_biom_string)
        self.mapping_file_headers = ['SampleID', 'HOST_SUBJECT_ID', 'Description']
        self.mapping_file_data = [['S1', 'A', 'candler fine sand'],
                                ['S2', 'A', 'candler fine sand'],
                                ['S3', 'A', 'candler fine sand'],
                                ['S4', 'B', 'candler fine sand'],
                                ['S5', 'B', 'candler fine sand'],
                                ['S6', 'B', 'candler fine sand'],
                                ['S7', 'C', 'candler fine sand'],
                                ['S8', 'C', 'candler fine sand']]
        self.expected_biom_string = output_biom_string_a

    def test_select_samples(self):

        # optimal case
        out_chosen_samples, out_biom_table = select_samples(\
            self.mapping_file_data, self.mapping_file_headers, self.biom_object,\
            40, 'HOST_SUBJECT_ID', 1, 1)

        self.assertEqual(['S7'], out_chosen_samples)
        self.assertEqual(out_biom_table,parse_biom_table(self.expected_biom_string))

        # undesired case
        with self.assertRaises(TableException):
            out_chosen_samples, out_biom_table = select_samples(\
                self.mapping_file_data, self.mapping_file_headers, 
                self.biom_object, 200, 'HOST_SUBJECT_ID', 1, 1)

input_biom_string = '{"rows": [{"id": "1", "metadata": null}, {"id": "2", "metadata": null}, {"id": "3", "metadata": null}, {"id": "4", "metadata": null}, {"id": "5", "metadata": null}, {"id": "6", "metadata": null}, {"id": "7", "metadata": null}, {"id": "8", "metadata": null}, {"id": "9", "metadata": null}, {"id": "10", "metadata": null}], "format": "Biological Observation Matrix 1.0.0", "data": [[0, 6, 14.0], [1, 3, 1.0], [1, 7, 5.0], [2, 4, 1.0], [2, 6, 4.0], [2, 7, 4.0], [4, 0, 4.0], [4, 1, 12.0], [4, 2, 4.0], [4, 3, 1.0], [4, 4, 1.0], [4, 5, 8.0], [4, 6, 8.0], [4, 7, 1.0], [5, 1, 11.0], [5, 2, 2.0], [5, 3, 2.0], [5, 4, 1.0], [5, 5, 4.0], [5, 6, 8.0], [5, 7, 22.0], [6, 1, 3.0], [6, 4, 2.0], [6, 5, 3.0], [6, 6, 4.0], [6, 7, 4.0], [7, 0, 2.0], [8, 0, 5.0], [8, 4, 1.0], [8, 6, 4.0], [8, 7, 1.0], [9, 6, 1.0]], "columns": [{"id": "S1", "metadata": null}, {"id": "S2", "metadata": null}, {"id": "S3", "metadata": null}, {"id": "S4", "metadata": null}, {"id": "S5", "metadata": null}, {"id": "S6", "metadata": null}, {"id": "S7", "metadata": null}, {"id": "S8", "metadata": null}], "generated_by": "BIOM-Format 1.0.0-dev", "matrix_type": "sparse", "shape": [10, 8], "format_url": "http://biom-format.org", "date": "2012-09-26T11:10:15.531807", "type": "OTU table", "id": null, "matrix_element_type": "float"}'
output_biom_string_a = '{"rows": [{"id": "1", "metadata": null}, {"id": "2", "metadata": null}, {"id": "3", "metadata": null}, {"id": "4", "metadata": null}, {"id": "5", "metadata": null}, {"id": "6", "metadata": null}, {"id": "7", "metadata": null}, {"id": "8", "metadata": null}, {"id": "9", "metadata": null}, {"id": "10", "metadata": null}], "format": "Biological Observation Matrix 1.0.0", "data": [[0, 0, 14.0], [2, 0, 4.0], [4, 0, 8.0], [5, 0, 8.0], [6, 0, 4.0], [8, 0, 4.0], [9, 0, 1.0]], "columns": [{"id": "S7", "metadata": null}], "generated_by": "test_biom", "matrix_type": "sparse", "shape": [10, 1], "format_url": "http://biom-format.org", "date": "2012-09-27T15:52:08.334245", "type": "OTU table", "id": null, "matrix_element_type": "float"}'

if __name__ == "__main__":
    main()
