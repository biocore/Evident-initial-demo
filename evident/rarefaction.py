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
from biom.table import TableException
from biom.parse import parse_biom_table
from qiime.colors import process_colorby
from qiime.collate_alpha import make_output_row
from qiime.make_rarefaction_plots import make_averages
from qiime.rarefaction import RarefactionMaker, get_rare_data
from qiime.util import compute_seqs_per_library_stats, FunctionWithParams
from qiime.parse import (parse_matrix, parse_mapping_file, parse_rarefaction,\
						parse_mapping_file_to_dict, parse_newick)

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
	"""calculates alpha rarefaction values but does not write output.
	this is a fileless implementation of single_file_alpha for use with 
	evident."""

	metrics_list = metrics.split(',')
	calcs = []
	for metric in metrics_list:
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
	# max_depth = num_seqs_per_sam due to earlier filtering
	# num_reps = num_iters
	rarefaction_step_size = int((max_depth - min_depth)/num_rare_depths)
	rare_maker = RarefactionMaker(biom_object, min_depth,  max_depth, 
		rarefaction_step_size, num_reps)
	rarefactions = rare_maker.rarefy_to_list() 
	return rarefactions

def generate_alpha_rarefaction_plots_from_point_in_omega(chosen_samples, map_file_tuple,
	filtered_biom_table, metrics_list, category, num_seqs_per_sam, num_iters, 
	tree_object=None, std_type='stddev'):
	"""Takes mapping file and biom table and generates alpha_rarefaction plots
	
	NOTES:
	Omega is the parameter space. It is N/{0}**3 where the axes are 
	number_of_subjects, number_of_samples_per_subject, and num_seqs_per_sam. 
	
	a nested dict with sampleIds of interest as keys, with values 
		equal to a list of as many dictionaries as there 
		number_of_axes. NOTE: result[sampleId1][0] gives the dictionary 
		for the avg, min, max of first axis of the pcoa, i.e. the 
		indexing is off by 1. 
		{'SampleId1:[{'avg': axis1, 'min':axis1, 'max':axis1}, 
					 {'avg': axis2, 'min':axis2, 'max':axis2}]
		 'SampleId2:[{'avg': axis1, 'min':axis1, 'max':axis1}...}
	outputs:
	"""
	# The minimum depth is defined by the size of the maximum depth
	num_reps = 4
	min_depth = int(ceil(num_seqs_per_sam / num_reps))

	rarefied_bioms = get_rarefactions(filtered_biom_table, min_depth, 
		num_seqs_per_sam, num_iters, num_reps)

	#convert metrics list to a string
	metrics = ','.join(metrics_list)
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
		color_prefs(map_file_tuple)

	# convert rarefaction strings/lists to proper format for make_average to use
	rares = {}
	for metric in all_metrics:
		ds = metrics_data_strings[metric]
		ds.append(col_header)
		rares[metric] = parse_rarefaction(ds)
		
	html_output = make_averages(prefs, data, background_color, label_color,\
		rares, 'dummy_fp', 75, 'png', None, False, std_type, output_type='memory')

	return html_output