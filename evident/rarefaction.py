#!/usr/bin/env python
# File created on 12 Jun 2012
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2012, Evident"
__credits__ = ["William Van Treuren", "Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"

# prevents runtime error with matplotlib, hence this statement must preceed
# any interaction with that module to avoid unwanted results
import os,tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

from math import ceil
from qiime.alpha_diversity import *
from qiime.colors import process_colorby
from qiime.rarefaction import RarefactionMaker
from qiime.collate_alpha import make_output_row
from qiime.parse import parse_matrix, parse_rarefaction

import logging

def build_color_preferences(mapping_file_tuple):
	"""defines and calculates color prefrences for alpha_rarefaction plots

	Inputs:
	mapping_file_tuple: parsed mapping file information i. e. (data, headers)

	Outputs:
	A list containing color_prefs, data, background_color, label_color for usage
	with make_averages or any other qiime visualization module
	"""
	data = {}

	# format the mapping file for use with process_by_color
	mapping_file_data, mapping_file_headers = mapping_file_tuple
	formatted_mapping_file =[mapping_file_headers]
	formatted_mapping_file.extend(mapping_file_data)
	data['map']= formatted_mapping_file

	background_color='black'
	label_color='white'

	# the colors will only be based on the data
	color_prefs, data = process_colorby(None, data, None)

	return color_prefs, data, background_color, label_color

def single_object_alpha(biom_object, metrics, tree_object):
	"""given a metric calculates alpha diversity of a biom object

	Inputs:
	biom_object: biom formatted OTU table
	metrics: list of alpha diversity metrics
	tree_object: tree object for the phylogenetic metrics

	Output:
	calculations: tab delimitted string with the calculations for the object
	"""

	calcs = []
	for metric in metrics:
		try:
			metric_f = get_nonphylogenetic_metric(metric)
			is_phylogenetic = False
		except AttributeError:
			try:
				metric_f = get_phylogenetic_metric(metric)
				is_phylogenetic = True
			except AttributeError:
				print 'something is wrong with conf file'
		c = AlphaDiversityCalc(metric_f, is_phylogenetic)
		calcs.append(c)

	all_calcs = AlphaDiversityCalcs(calcs)

	result = all_calcs(data_path=biom_object, tree_path=tree_object,
		log_path=None)
	return all_calcs.formatResult(result)

def get_rarefactions(biom_object, minimum, maximum, iterations, steps):
	"""rarify biom object and return rarefactions

	Inputs:
	biom_object: object to rarefy
	minimum: starting point for the rarefactions
	maximum: ending point for the rarefactions
	iterations: repetitions per rarefaction depth
	steps: number of levels between minimum and maximum

	Outputs:
	list of 3 element tuples, where each tuple contains as a 1st element the
	rarefaction depth, as a 2nd element the iteration number and as a 3rd
	element the rarefied biom corresponding to this depth
	"""

	rarefaction_step_size = int((maximum - minimum)/steps)

	rarefaction_maker = RarefactionMaker(biom_object, minimum,  maximum,\
		rarefaction_step_size, iterations)
	rarefactions = rarefaction_maker.rarefy_to_list()

	return rarefactions

def _format_rarefactions(rarefaction_data, samples):
	"""utility to get rarefaction data into a make_averages compatible format

	Inputs:
	rarefaction_data: dictionary of alpha-div metrics and their values
	samples: list of sample identifiers for each of the alpha-div calculations

	Output:
	output: qiime.make_rarefaction_plots.make_averages compatible rarefaction 
	data
	"""
	output = {}
	# each key is a metric, each value is a list of alpha diversity values
	for key, value in rarefaction_data.iteritems():

		# build the header of the table
		_buffer = [['', 'sequences per sample', 'iteration']]
		_buffer[0].extend(samples)
		_buffer.append([])

		# extract each of the rarefaction labels
		_buffer.append([row[0] for row in value])

		# all the rarefaction values must be formatted as floats
		full_matrix = []
		for row in value:
			matrix_row = []
			for element in row[1::]:
				# some elements are strings and some elements are floats
				try:
					matrix_row.append(float(element))
				except TypeError:
					matrix_row.append(element)
			full_matrix.append(matrix_row)
		_buffer.append(full_matrix)

		output[key] = tuple(_buffer)

	return output

def generate_alpha_rarefaction_data_from_point_in_omega(biom_object, metrics,
													sequences, iterations,
													tree_object=None):
	"""generate alpha rarefaction data from a biom table and mapping file

	Inputs:
	biom_object: OTU table to be rarefied and used to compute alpha diversity
	metrics: list of metrics, phylogenetic or non phylogenetic
	sequences: maximum number of sequences for the rarefaction plots
	iterations: number of repetitions per rarefaction
	tree_object: tree to perform the phylogenetic operations, default is None

	Output:
	alpha_rarefaction_data: dictionary where the keys are alpha diversity
	metrics and the values are tuples; in these tuples the first element is a
	list of column headers for an alpha diversity file, the second element is a
	list of row headers for an alpha diversity file and the third element is a
	list of lists containing the alpha diversity data computed at multiple
	rarefaction depths and as many iterations as specified.
	"""
	# The minimum depth is defined by the size of the maximum depth
	steps = 4
	min_depth = int(ceil(sequences / steps))

	# get a rarefied biom with the proper identifiers
	rarefied_bioms_list = get_rarefactions(biom_object, min_depth, sequences,\
		iterations, steps)

	alpha_rs = {}
	alpha_filenames = []
	# rarefy all the biom objects and get the alpha diversity values
	for rarefied_biom in rarefied_bioms_list:
		# this tag contains data about the iteration and the depth
		identifier = 'alpha_rare_%s_%s' % (str(rarefied_biom[0]), str(rarefied_biom[1]))
		alpha_values = single_object_alpha(rarefied_biom[2], metrics, tree_object)
		alpha_rs[identifier] = (rarefied_biom[0], rarefied_biom[1], alpha_values.split('\n'))
		alpha_filenames.append(identifier)

	# use the rarefaction with the fewest sequences per sample as the reference
	ref_rare = single_object_alpha(rarefied_bioms_list[0][2], metrics,\
		tree_object=tree_object).split('\n')
	all_metrics, all_samples, example_data = parse_matrix(ref_rare)

	# build a dictionary with the data for each of the metrics specified
	metrics_data = {}
	for metric in all_metrics:
		per_metric_data = []
		for filename in alpha_filenames:
			f_metrics, f_samples, f_data = parse_matrix(alpha_rs[filename][2])
			per_metric_data.append(make_output_row(f_metrics, metric,\
				f_samples, f_data, filename, len(all_samples), all_samples))
		metrics_data[metric] = per_metric_data

	# now format the dictionary to make it compatible with make_averages
	alpha_rarefaction_data = _format_rarefactions(metrics_data, all_samples)

	return alpha_rarefaction_data
