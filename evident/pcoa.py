#!/usr/bin/env python
# File created on 12 Jun 2012
#from __future__ import division

__author__ = "Antonio Gonzalez"
__copyright__ = "Copyright 2012, Evident"
__credits__ = ["Antonio Gonzalez, William Van Treuren"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "Antonio Gonzalez"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

"""code for doing pcoa analysis"""

from numpy import array
from cStringIO import StringIO

from qiime.rarefaction import get_rare_data
from qiime.principal_coordinates import pcoa
from qiime.beta_diversity import single_object_beta
from qiime.parse import parse_mapping_file, mapping_file_to_dict, parse_coords
from emperor.format import (format_pcoa_to_js, format_mapping_file_to_js, 
    format_taxa_to_js, format_vectors_to_js, format_emperor_html_footer_string, 
    format_comparison_bars_to_js, EMPEROR_HEADER_HTML_STRING)
from emperor.util import preprocess_coords_file


def make_pcoa_plot(pcoa_headers, pcoa_files, eigenvalues, coords_pct, map_headers, 
        map_data, coords_low=None, coords_high=None, jackkifing_controls=False): 
    """ generate pcoa plot, this is based on make_emperor.py 

    Input:
    pcoa_headers: list of headers from the pcoa file
    pcoa_files: numpy array with the values of the pcoa
    eigenvalues: list of the eigenvalues of the pcoa
    coords_pct: list of the values of the coords points
    map_headers: list of strings with the mapping headers
    map_data: list of lists with the values of the mapping file
    
    Output:
    WebGL string representing the PCoA plot
    """

    webgl_string = EMPEROR_HEADER_HTML_STRING
    webgl_string += format_mapping_file_to_js(map_data, map_headers, map_headers)
    webgl_string += format_pcoa_to_js(pcoa_headers, pcoa_files,
            eigenvalues, coords_pct, custom_axes=None, coords_low=coords_low,
            coords_high=coords_high, number_of_axes=10, number_of_segments=8)
    
    webgl_string += format_taxa_to_js([], [], [])
    webgl_string += format_vectors_to_js(map_data, pcoa_files, map_headers,
        pcoa_headers, None, None)
    webgl_string += format_comparison_bars_to_js(pcoa_files, pcoa_headers, 0, False)
    webgl_string += format_emperor_html_footer_string(False, jackkifing_controls, \
        False, False)
    
    return webgl_string

def generate_pcoa_cloud_from_point_in_omega(map_headers, map_data, biom_object, metric, 
        sequences, iterations, axes, tree_object=None):
    """run the randomisations and get a WebGL PCoA plot string representation

    Input:
    mapping_file_tuple: data and headers tuple for representing the mapping file
    biom_object: otu table biom object
    metric: string of the name for the beta diversity metric, i. e. 'unifrac'
    sequences: number of sequences per sample
    iterations: number of iterations to generate the pcoa plot
    axes: number of axes to account for
    tree_object: tree to perform the beta diversity calculation

    Output:
    WebGL string representing the PCoA plot
    """
    
    pcoa_input = {'pcoa_headers':[], 'pcoa_values':[], 'eigenvalues':[], 'coords_pct':[]}
    for i in range(iterations):
        rare_biom_table = get_rare_data(biom_object, sequences)
        beta_dm = single_object_beta(rare_biom_table, metric, tree_object)
        pcoa_results = pcoa(beta_dm)

        pcoa_file = StringIO()
        pcoa_file.write(pcoa_results)
        pcoa_file.seek(0)
        pcoa_headers, pcoa_values, eigenvalues, coords_pct = parse_coords(pcoa_file)
        pcoa_file.close()
        pcoa_input['pcoa_headers'].append(pcoa_headers)
        pcoa_input['pcoa_values'].append(pcoa_values)
        pcoa_input['eigenvalues'].append(eigenvalues)
        pcoa_input['coords_pct'].append(coords_pct)
    
    if iterations==1:
        coords_headers = pcoa_input['pcoa_headers'][0]
        coords_data = pcoa_input['pcoa_values'][0]
        coords_eigenvalues = pcoa_input['eigenvalues'][0]
        coords_pct = pcoa_input['coords_pct'][0]
        coords_low, coords_high = None, None
    else:
        coords_headers, coords_data, coords_eigenvalues, coords_pct, coords_low,\
            coords_high, clones = preprocess_coords_file(pcoa_input['pcoa_headers'],
            pcoa_input['pcoa_values'], pcoa_input['eigenvalues'], 
            pcoa_input['coords_pct'], map_headers, map_data, custom_axes=None, 
            jackknifing_method='IQR', is_comparison=False)
    
    return make_pcoa_plot(coords_headers, coords_data, coords_eigenvalues, coords_pct, \
        map_headers, map_data, coords_low, coords_high, True)