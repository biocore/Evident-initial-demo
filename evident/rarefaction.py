#!/usr/bin/env python
# File created on 12 Jun 2012
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2012, Evident"
__credits__ = ["William Van Treuren"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "William Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"

# prevents runtime error with matplotlib
import os,tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

from math import ceil
from qiime.alpha_diversity import *
from qiime.colors import process_colorby
from qiime.rarefaction import RarefactionMaker
from qiime.collate_alpha import make_output_row
from qiime.make_rarefaction_plots import make_averages
from qiime.parse import parse_matrix, parse_rarefaction

import logging

def color_prefs(parsed_mf):
	"""defines and calculates color prefrences for alpha_rarefaction plots"""
	data = {}
	mapping,headers = parsed_mf # mapping file is tuple, doesnt have comments
	new_mapping=[headers]
	new_mapping.extend(mapping)
	data['map']=new_mapping
	background_color='black'
	label_color='white'
	ball_scale=1.0
	arrow_colors={'line_color':'white', 'head_color':'red'}
	color_prefs, data = process_colorby(colorby=None, data=data,
		color_prefs=None)
	return color_prefs, data, background_color, label_color, ball_scale, \
		arrow_colors

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

def get_rarefactions(biom_object,min_depth,max_depth,num_reps,num_rare_depths):
	"""rarify biom object and return rarefactions"""
	rarefaction_step_size = int((max_depth - min_depth)/num_rare_depths)
	rare_maker = RarefactionMaker(biom_object, min_depth,  max_depth, 
		rarefaction_step_size, num_reps)
	rarefactions = rare_maker.rarefy_to_list() 
	return rarefactions

def _format_rarefactions(rarefaction_data, samples):
	""""""


def generate_alpha_rarefaction_plots_from_point_in_omega(mapping_file_tuple,
														biom_object, metrics, 
														sequences, iterations, 
														tree_object=None, 
														std_type='stddev'):
	"""generate alpha rarefaction plots from a biom table and mapping file

	Inputs:
	mapping_file_tuple: mapping file data and headers in a tuple (data, headers)
	biom_object: OTU table in a biom format corresponding to mapping_file_tuple
	metrics: list of metrics, phylogenetic or non phylogenetic
	sequences: maximum number of sequences for the rarefaction plots
	iterations: number of repetitions per rarefaction
	tree_object: tree to perform the phylogenetic operations, default is None
	std_type: calculation to perform for the error bars, can be standard
	deviation (stddev) or standard error (stderr), default is stddev 

	Outputs:
	html_string: HTML formatted string with the rarefaction plots for the given
	parameters
	"""
	# The minimum depth is defined by the size of the maximum depth
	steps = 4
	min_depth = int(ceil(sequences / steps))

	rarefied_bioms = get_rarefactions(biom_object, min_depth, 
		sequences, iterations, steps)

	alpha_rs = {}
	alpha_filenames = []
	for rb in rarefied_bioms:
		key = 'alpha_rare_%s_%s' % (str(rb[0]), str(rb[1]))
		alpha_values = single_object_alpha(rb[2], metrics, tree_object)
		alpha_rs[key] = (rb[0], rb[1], alpha_values.split('\n'))
		alpha_filenames.append(key)

	# use the rarefaction with the fewest seqs/sample as the refrence 
	ref_rare = single_object_alpha(rarefied_bioms[0][2], metrics, 
		tree_object=tree_object).split('\n')
	all_metrics, all_samples, example_data = parse_matrix(ref_rare)

	num_cols = len(all_samples)

	#REFACTOR
	# metrics_data is dict:
	# {'Shannon':['alpha_rare_seqs_20_iter_0','5.64545',...], 'Chao1':...} 
	metrics_data = {}
	for metric in all_metrics:
		metric_file_data = []
		for fname in alpha_filenames:
			f_metrics, f_samples, f_data = parse_matrix(alpha_rs[fname][2])
			metric_file_data.append(\
			  make_output_row(f_metrics, metric, f_samples, 
				f_data, fname,num_cols,all_samples))
		metrics_data[metric] = metric_file_data

	# convert the metrics_data to one long string to fool parse_rarefaction.
	# need to be careful; parse_rarefaction thinks that its the column header
	# if it starts with a tab.	make column header start with a tab, and 
	# everything else start without a tab. include sequences and iteration data

	# REFACTOR
	metrics_data_strings = {}
	for metric in metrics_data:
		metrics_lines = []
		for line in metrics_data[metric]:
			# need the following form alpha_rare_seqs_20_iter_0\t20\t0
			list_of_components = line[0].split('_')
			iter = list_of_components.pop()
			num_seqs = list_of_components.pop()
			stringified_line = line[0]+'\t'#+num_seqs+'\t'+iter+'\t'

			for list_element in line[1:]: # line[0] is the name, used above
				stringified_line = stringified_line + str(list_element) + '\t'
			# get rid of the final tab
			stringified_line = stringified_line[:-1]
			metrics_lines.append(stringified_line)
		metrics_data_strings[metric] = metrics_lines

	# now we need to create the column header line	  
	all_samples_string = '\t'.join(all_samples)
	col_header = '\tsequences per sample\titeration\t'+all_samples_string +'\n'
	# need the \n to match from file

	prefs, data, background_color, label_color, ball_scale, arrow_colors = \
		color_prefs(mapping_file_tuple)

	# convert rarefaction strings/lists to proper format for make_average to use
	rares = {}
	for metric in all_metrics:
		ds = metrics_data_strings[metric]
		ds.append(col_header)
		rares[metric] = parse_rarefaction(ds)
		
	html_output = make_averages(prefs, data, background_color, label_color,\
		rares, 'dummy_fp', 75, 'png', None, False, std_type, output_type='memory')

	return html_output