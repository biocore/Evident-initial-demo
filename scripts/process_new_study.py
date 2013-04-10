#!/usr/bin/env python
# File created on 01 Apr 2013
# This script is based on: beta_diversity_through_plots.py, alpha_rarefaction.py
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2012, Evident"
__credits__ = ["Antonio Gonzalez Pena", "Will Van Treuren", "Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"


from biom.parse import parse_biom_table
from qiime.filter import filter_mapping_file, filter_samples_from_otu_table
from qiime.util import parse_command_line_parameters, make_option, get_options_lookup
from qiime.parse import parse_qiime_parameters, parse_mapping_file
from qiime.workflow.util import (print_commands,
                                 call_commands_serially,
                                 print_to_stdout,
                                 no_status_updates,
                                 validate_and_set_jobs_to_start)
from qiime.workflow.downstream import run_beta_diversity_through_plots, run_alpha_rarefaction
from qiime.util import load_qiime_config, parse_command_line_parameters, get_options_lookup
from qiime.format import format_biom_table, format_mapping_file
from evident.map_sample_space import get_sorted_counts_per_sample, make_selectors
from os import makedirs
from numpy import inf
from shutil import copyfile

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Processes a biom and mapping file downloaded from " +\
    "the QIIME database to be added to Evident"
script_info['script_description'] = """This script takes a biom and mapping file \
downloaded from the database, together with how many sequences per sample and the \
subject id in the mapping file and processes the study to be added to Evident."""
script_info['script_usage'] = [("Process a new study","""To process a new study you need \
to download a study from the database, biom and mapping file, select a rarefaction level \
to perform the analyses, and define the column in the mapping file that has the unique \
identifier of the subjects (for example: HOST_SUBJECT_ID) and then run the following \
command:""", """%prog -i otu_table.biom -m mapping_file.txt -o processed_study -e 1000 \
-s HOST_SUBJECTY""")]
script_info['script_usage'].append(("Process a new study in parallel","""To process a \
study in parallel, using 10 jobs, you can use this command:""", """%prog -i \
otu_table.biom -m mapping_file.txt -o processed_study -e 1000 -s HOST_SUBJECTY -aO 10"""))
script_info['output_description']="""The script creates a raw.biom (original file), \
an even sampled biom file, a selectors.txt file that contains information of how evident \
should behave in the main GUI, a cleaned mapping file, a study_preference file that has \
some basic information about the study, and alpha & beta calculations. """
script_info['required_options'] = [\
 make_option('-i','--otu_table_fp',type='existing_filepath',
            help='the input biom table [REQUIRED]'),
 make_option('-m','--mapping_fp',type='existing_filepath',
            help='path to the mapping file [REQUIRED]'),
 make_option('-o','--output_dir',type='new_dirpath',
            help='the output directory [REQUIRED]'),
 make_option('-e','--seqs_per_sample',type="int",help='the input filepath'),\
 make_option('-s','--subject_name',help='the name of the subject category in the ' + 
            'mapping file, i.e. "INDIVIDUAL", "HOST_INDIVIDUAL"'),\
]
script_info['optional_options'] = [\
 make_option('-t','--tree_fp',type='existing_filepath',
        help='path to the tree file [default: %default]',
        default="/evident/data/gg_97_otus_4feb2011.tre"),
 make_option('-p','--parameter_fp',type='existing_filepath',
    help='path to the parameter file, which specifies changes'+\
        ' to the default behavior. '+\
        'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .'+\
        ' [if omitted, default values will be used]'),
 make_option('-f','--force',action='store_true',
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),
 make_option('-w','--print_only',action='store_true',
        dest='print_only',help='Print the commands but don\'t call them -- '+\
        'useful for debugging [default: %default]',default=False),
 make_option('-a','--parallel',action='store_true',
        dest='parallel',default=False,
        help='Run in parallel where available [default: %default]'),
 options_lookup['jobs_to_start_workflow'],
 make_option('-n','--min_seqs_sample',type=int,help='min number of sequence to ' + 
        'be considered a valid state, default=%default',default=100),\
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    otu_table_fp = opts.otu_table_fp
    output_dir = opts.output_dir
    mapping_fp = opts.mapping_fp
    tree_fp = opts.tree_fp
    verbose = opts.verbose
    print_only = opts.print_only
    seqs_per_sample = int(opts.seqs_per_sample)
    parallel = opts.parallel
    min_seqs_sample = opts.min_seqs_sample
    subject_category = opts.subject_name

    try:
        makedirs(output_dir)
    except OSError:
        if opts.force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            option_parser.error("Output directory already exists. Please choose"
                " a different directory, or force overwrite with -f.")


    ## ******************** make_evident_selectors ********************
    ## The code for make_evident_selectors.py is here and has to go before the params
    ## validation as we need to know the main cats before creating the params file
    map_data, headers, comments = parse_mapping_file(open(mapping_fp, 'U'))
    biom_table = parse_biom_table(open(otu_table_fp, 'U'))

    # getting valid samples from biom file
    real_map_headers, real_map_data = filter_mapping_file(map_data, headers,\
        biom_table.SampleIds, include_repeat_cols=False)

    if subject_category not in real_map_headers:
        raise ValueError, 'This column: %s is not in the mapping file, try %s'%\
            (subject_category, real_map_headers)

    sorted_counts_per_sample = get_sorted_counts_per_sample(biom_table)

    mapping_file_tuple = (real_map_data, real_map_headers)

    # calculate the available subjects at each rarefaction level
    results, main_map_cat = make_selectors(sorted_counts_per_sample, min_seqs_sample,\
        mapping_file_tuple, subject_category, verbose=verbose)

    fout = open(output_dir + '/selectors.txt','w')
    fout.write('#Sequences\tSubjects\tSamples\tMetadata\n')
    fout.write('\n'.join(results))
    fout.close()
    
    fout = open(output_dir + '/mapping_file.txt','w')
    fout.write(format_mapping_file(real_map_headers, real_map_data))
    fout.close()
    ## ******************** make_evident_selectors ********************

    fout = open(output_dir + '/study_preferences.txt','w')
    fout.write('%d\n' % seqs_per_sample)
    fout.write('%s\n' % subject_category)
    fout.close()

    ## ******************** filter_samples_from_otu_table ********************
    ## Filtering original biom file to only have samples above the max length to avoid
    ## ugly plots
    alpha_biom_file = output_dir + '/filtered_otu_table_for_alpha.biom'
    fout = open(alpha_biom_file,'w')
    sample_ids_to_keep = biom_table.SampleIds
    filtered_otu_table = filter_samples_from_otu_table(biom_table,
                                                       sample_ids_to_keep,
                                                       min_count=seqs_per_sample,
                                                       max_count=inf)
    fout.write(format_biom_table(filtered_otu_table))
    fout.close()
    ## ******************** filter_samples_from_otu_table ********************

    if opts.parameter_fp:
        try:
            parameter_f = open(opts.parameter_fp, 'U')
        except IOError:
            raise IOError,\
             "Can't open parameters file (%s). Does it exist? Do you have read access?"\
             % opts.parameter_fp
        params = parse_qiime_parameters(parameter_f)
        parameter_f.close()
    else:
        params = parse_qiime_parameters(
            ['beta_diversity:metrics unweighted_unifrac',\
             'make_rarefaction_plots:prefs_path %s' % output_dir + '/prefs.txt',
             'make_rarefaction_plots:colorby %s' % ','.join(main_map_cat), 
             'make_rarefaction_plots:output_type memory', 
             'multiple_rarefactions:min %d' % int(seqs_per_sample/4),
             'multiple_rarefactions:max %d' % (seqs_per_sample+1),
             'multiple_rarefactions:step %d' % int(seqs_per_sample/4),
             'multiple_rarefactions:num-reps 4',
            ])
        # empty list returns empty defaultdict for now
    
    jobs_to_start = opts.jobs_to_start
    default_jobs_to_start = qiime_config['jobs_to_start']
    validate_and_set_jobs_to_start(params,
                                   jobs_to_start,
                                   default_jobs_to_start,
                                   parallel,
                                   option_parser)


    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially
    
    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates
    
    copyfile(otu_table_fp, output_dir + '/raw.biom')
    
    run_beta_diversity_through_plots(otu_table_fp=otu_table_fp,
     mapping_fp=mapping_fp,
     output_dir=output_dir,
     command_handler=command_handler,
     params=params,
     qiime_config=qiime_config,
     color_by_interesting_fields_only=False,
     sampling_depth=seqs_per_sample,
     histogram_categories=None,
     tree_fp=tree_fp,
     parallel=parallel,
     suppress_3d_plots=True,
     suppress_2d_plots=True,
     status_update_callback=status_update_callback)
    
    output_dir = output_dir + '/alpha'
    run_alpha_rarefaction(otu_table_fp=alpha_biom_file,\
     mapping_fp=mapping_fp,\
     output_dir=output_dir,\
     command_handler=command_handler,\
     params=params,
     qiime_config=qiime_config,\
     tree_fp=tree_fp,\
     num_steps=4,\
     parallel=parallel,\
     min_rare_depth=10,
     max_rare_depth=20,
     status_update_callback=status_update_callback,
     plot_stderr_and_stddev=True)
    
     

if __name__ == "__main__":
    main()