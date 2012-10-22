#!/usr/bin/env python
# File created on 07 Mar 2012
from __future__ import division

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2012, Evident"
__credits__ = ["Will Van Treuren", "Antonio Gonzalez Pena",
"Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"


from biom.parse import parse_biom_table
from qiime.parse import parse_mapping_file
from qiime.filter import filter_mapping_file
from qiime.util import parse_command_line_parameters, make_option
from evident.map_sample_space import (get_sorted_counts_per_sample, make_selectors)

script_info = {}
script_info['brief_description'] = "Create a file with the ranges a user can select using the Evident GUI"
script_info['script_description'] = "This script looks at the available "+\
"samples of a study and the number of sequences for each sample, then creates"+\
" a tab delimited file with four columns: number of sequences (#Sequences), "+\
"number of subjects (Subjects), number of samples (Samples) and the " +\
"metadata available (Metadata). To contabilize the subjects, you will need "+\
"choose with the unique-per-subject-identifier."
script_info['script_usage'] = [("Creation of a selectors file for Evident","To be able to use a dataset with Evident, you will need to have an BIOM formatted OTU table, a QIIME formatted mapping file and an Evident formatted selectors file created by this script: ","script_map_sample_space.py -m mapping.txt -i otu_table.biom -s HOST_SUBJECT_ID  -o evident_selectors.txt")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-m','--mapping',type="existing_filepath",help='the input mapping file'),\
 make_option('-i','--biom_file',type="existing_filepath",help='the input biom file'),\
 make_option('-s','--subject_name',help='the name of the subject category in the mapping file, i.e. "INDIVIDUAL", "HOST_INDIVIDUAL"'),\
 make_option('-o','--clean_fp',type="new_filepath",help='the name of the clean file path'),\
]
script_info['optional_options'] = [\
# Example optional option
 make_option('-n','--min_seqs_sample',type=int,help='min number of sequence to be considered '+\
            'a valid state, default=%default',default=100),\
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    map_fp = opts.mapping
    biom_fp = opts.biom_file
    min_seqs_sample = opts.min_seqs_sample
    subject_category = opts.subject_name

    cleaned_fp = opts.clean_fp
    verbose = opts.verbose

    map_data, headers, comments = parse_mapping_file(open(map_fp, 'U'))
    biom_table = parse_biom_table(open(biom_fp, 'U'))

    # getting valid samples from biom file
    real_map_headers, real_map_data = filter_mapping_file(map_data, headers,\
        biom_table.SampleIds, include_repeat_cols=False)

    if subject_category not in real_map_headers:
        raise ValueError, 'This column: %s is not in the mapping file, try %s'%\
            (subject_category, real_map_headers)

    sorted_counts_per_sample = get_sorted_counts_per_sample(biom_table)

    mapping_file_tuple = (real_map_data, real_map_headers)

    # calculate the available subjects at each rarefaction level
    results = make_selectors(sorted_counts_per_sample, min_seqs_sample,\
        mapping_file_tuple, subject_category, verbose=verbose)

    # save the output
    fout = open(cleaned_fp,'w')
    fout.write('#Sequences\tSubjects\tSamples\tMetadata\n')
    fout.write('\n'.join(results))
    fout.close()

if __name__ == "__main__":
    main()
