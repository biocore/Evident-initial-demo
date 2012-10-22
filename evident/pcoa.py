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
from qiime.parse import parse_mapping_file, mapping_file_to_dict

def make_pcoa_webgl_string_from_files(pcoa_headers, points, map_headers, map_data):
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

    webgl_string = make_pcoa_plot(ellipses, (map_data, map_headers))
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

    Assume that the following represents the x pcoa axis, and each + is the
    x coordinate of a point produced by the pcoa plots code.

    x_axis  = +...+..+.+.+.+..........0..........+.+.+..+..............+
    density =     ---------------------------------------
    Now, in the ideal world, we would want the ellipse to be darker, or the
    transparency to change based on the fact that even though we might every
    once in a while get an outlier point which would be either really high or
    really low on the x axis, the majority of the probability density is located
    in a certain interval. thus, when we need to calculate the length of the
    ellipsoid axes we have several possible ways we could go: 
        1. Take the average of the abs of the x points and then that is the 
        ellipse x axis. This means that the ellipse we would produce would not
        actually show the full extent of how far the data could move, but would
        weight the ellipse.
        2. Take a 95 percent confidence interval, pick a length for the axis
        that includes at least 95 percent of the x points.
        3. Weighting function
        4. Use the max CURRENTLY USING 1
    """
    sampleId_to_coords = {}
    sampleId_center_and_axes = {}
    for pcoa_string in sampled_pcoa_strings:
        for line in pcoa_string.split('\n'):
            line = line.split('\t')
            if line[0] in sampleIds:
                if line[0] not in sampleId_to_coords.keys():
                    sampleId_to_coords[line[0]] = {'coords':[[] for i in range(1,number_of_axes+1,1)]}
                    sampleId_center_and_axes[line[0]] = { 'center':[], 'axes_radii':[] }
                for axis in range(number_of_axes):
                    sampleId_to_coords[line[0]]['coords'][axis].append(float(line[axis+1]))

    for samId in sampleId_center_and_axes:
        for axis in range(number_of_axes):
            center = array(sampleId_to_coords[samId]['coords'][axis]).mean()
            sampleId_center_and_axes[samId]['center'].append(center)

            dfc = abs(\
                array(sampleId_to_coords[samId]['coords'][axis]) - \
                      sampleId_center_and_axes[samId]['center'][axis]).mean()
            sampleId_center_and_axes[samId]['axes_radii'].append(dfc)

    return sampleId_center_and_axes

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
    ellipse_coords_by_sampleId = get_pcoa_ellipsoid_coords(pcoa_list, axes,\
        full_id_list)

    # check the ellipses are created correctly
    if type(ellipse_coords_by_sampleId) == type(''):
        raise ValueError, 'Could not create PCoA plot'

    webgl_string = make_pcoa_plot(ellipse_coords_by_sampleId, mapping_file_tuple)
    return webgl_string
