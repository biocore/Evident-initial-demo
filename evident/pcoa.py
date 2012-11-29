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

from evident import make_pcoa_plot

from qiime.rarefaction import get_rare_data
from qiime.principal_coordinates import pcoa
from qiime.beta_diversity import single_object_beta
from qiime.parse import parse_mapping_file, mapping_file_to_dict, parse_coords

def make_pcoa_webgl_string_from_files(pcoa_headers, points, pct_var, map_headers, map_data):
    """make a WebGL representation of a pcoa plot with the provided data

    Input:
    pcoa_headers: list of headers of a pcoa file
    points: coordinates of a pcoa file
    map_headers: list of headers for a mapping file
    map_data: list of lists containing the data of a mapping file

    Output:
    webgl_string: string representation of the pcoa plot for WebGL display
    """
    ellipses = {}
    for h,p in zip(pcoa_headers, points):
        ellipses[h] = { 'axes_radii': [0.0, 0.0, 0.0], 'center': p[:3] }

    # check the ellipses are created correctly
    if ellipses == {}:
        raise ValueError, 'Could not create PCoA plot'

    webgl_string = make_pcoa_plot(ellipses, (map_data, map_headers), pct_var[:3])
    return webgl_string

def get_pcoa_ellipsoid_coords(sampled_pcoa_strings, number_of_axes, sampleIds):
    """gets min, max, average for each axis in the number_of_axes

    Inputs:
    sampled_pcoa_strings: list of PCoA strings
    number_of_axes: number of axes for the PCoA
    SampleIDs: list of the contained identifiers

    Outputs:
    A nested dict with sampleIds of interest as keys, with values equal to a
    list of as many dictionaries as there number_of_axes.

    NOTE: result[sampleId1][0] gives the dictionary for the avg, min, max of
    first axis of the pcoa, i.e. the indexing is off by 1.
    {'SampleId1:{'center':[x,y,z,w]
                 'axes_radii':[x_radius, y_radius, z_radius, w_radius]},
     'SampleId2:{'center':[x,y,z,w]
                 'axes_radii':[x_radius, y_radius, z_radius, w_radius]}

    The ellipsoids that are calculated by this function are NOT minimum spanning
    ellipsoids. This function looks at all the iterations for a given sample
    and the calculates the center and max distance from the center of each of 
    those iterations along each individual axis. While this ensures that the 
    most extreme points on any given axis are contained within the ellipse it 
    isn't guaranteed that every point will be contained. As an example, 
    circumscribe a rectangle about an ellipse so that the ellipse and rectangle 
    are co-linear at 4 points (i.e. rectangle is just big enough to contain 
    the ellipse). Take a corner of the rectangle and move it towards the center
    of the ellipse so that its x and y components diminish in magnitude. This 
    new point will not be contained within the ellipse, but its x and y 
    coordinates will be less extreme than the most extreme points of the 
    ellipse (the points that are colinear with the minor and major axis.)
    """
    sampleId_to_coords = {}
    sampleId_center_and_axes = {}
    for pcoa_string in sampled_pcoa_strings:
        coord_header, coords, eigvals, pct_var = parse_coords(pcoa_string.split('\n'))
        if 'variation explained' not in sampleId_to_coords:
            sampleId_to_coords['variation explained'] = []
        sampleId_to_coords['variation explained'].append(pct_var[:number_of_axes])
        
        for sampleName, values in zip(coord_header, coords):
            if sampleName not in sampleId_to_coords.keys():
                sampleId_to_coords[sampleName] = {'coords':[[] for i in range(1,number_of_axes+1,1)]}
                sampleId_center_and_axes[sampleName] = { 'center':[], 'axes_radii':[] }
            for axis in range(number_of_axes):
                sampleId_to_coords[sampleName]['coords'][axis].append(values[axis])
    sampleId_to_coords['variation explained'] = array(sampleId_to_coords['variation explained']).mean(0)
    
    for samId in sampleId_center_and_axes:
        for axis in range(number_of_axes):
            center = array(sampleId_to_coords[samId]['coords'][axis]).mean()
            sampleId_center_and_axes[samId]['center'].append(center)

            dfc = abs(\
                array(sampleId_to_coords[samId]['coords'][axis]) - \
                      sampleId_center_and_axes[samId]['center'][axis]).mean()
            sampleId_center_and_axes[samId]['axes_radii'].append(dfc)

    return sampleId_center_and_axes, sampleId_to_coords

def generate_pcoa_cloud_from_point_in_omega(mapping_file_tuple, biom_object,
                                            metric, sequences, iterations, axes,
                                            tree_object=None):
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

    # get a list of the SampleIds
    full_id_list = mapping_file_to_dict(mapping_file_tuple[0], mapping_file_tuple[1]).keys()

    pcoa_list = []
    for i in range(iterations):
        rare_biom_table = get_rare_data(biom_object, sequences)
        beta_dm = single_object_beta(rare_biom_table, metric, tree_object)
        pcoa_results = pcoa(beta_dm)

        pcoa_list.append(pcoa_results)

    # convert the list of pcoa lines into ellipsoid coords
    ellipse_coords_by_sampleId, sampleId_to_coords  = get_pcoa_ellipsoid_coords(pcoa_list, axes, full_id_list)
        
    # check the ellipses are created correctly
    if type(ellipse_coords_by_sampleId) == type(''):
        raise ValueError, 'Could not create PCoA plot'

    webgl_string = make_pcoa_plot(ellipse_coords_by_sampleId, mapping_file_tuple, sampleId_to_coords['variation explained'])
    return webgl_string
