#!/usr/bin/env python
# File created on 07 Mar 2012
from __future__ import division

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2012, Evident"
__credits__ = ["Will Van Treuren", "Antonio Gonzalez Pena", "Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"

from random import shuffle
from biom.exception import TableException
from qiime.rarefaction import get_rare_data

import logging

def select_samples(map_data, headers, biom_table, depth, unique_id_column, 
                    subjects, samples_per_subject):
    """
    Randomly select a list of IDs with enough sequeneces per subject

    Input:
    map_data: rows of the mapping file without the headers and the comments
    headers: headers of the mapping file
    biom_table: table object 
    depth: number of sequences pers sample
    unique_id_column: column header to identify unique subjects in the mapping 
    file i.e. HOST_SUBJECT_ID
    subjects: number of subjects to include in the resulting samples
    samples_per_subject: number of samples pers subject to include in the
    resulting samples

    Output:
    chosen_samples: a list of SampleIds with as many lists as subjects where 
    each list has as many elmenents as samples per subject
    final_biom_table: a biom table object containing only the 'chosen_samples'
    """

    unique_id_column_index = headers.index(unique_id_column)
    rare_biom_table = get_rare_data(biom_table, depth)

    # make a dictionary of each subject with its corresponding list of SampleIds
    per_subject_sample_ids = {}
    for row in map_data:
        if row[0] not in rare_biom_table.SampleIds:
            continue

        if  row[unique_id_column_index] not in per_subject_sample_ids:
            per_subject_sample_ids[row[unique_id_column_index]] = []
        per_subject_sample_ids[row[unique_id_column_index]].append(row[0])

    # subsampling samples per individual
    subsampled_ids = {}
    subject_keys = []
    for k,v in per_subject_sample_ids.items():
        if len(v)>=samples_per_subject:
            subsampled_ids[k] = v[:samples_per_subject]
            subject_keys.append(k)

    ##open('/tmp/tmp.txt','a').write('%s %d %s\n\n' % (k, len(subsampled_ids[k]), subsampled_ids[k]))

    # subsampling subjects
    shuffle(subject_keys)
    chosen_samples = []
    for k in subject_keys[:subjects]:
        chosen_samples.extend(subsampled_ids[k])
    
    # creating new biom file with only the good samples
    try:
        final_biom_table = biom_table.filterSamples(lambda v,id,md: id in chosen_samples)
    except TableException:
        raise TableException, "Using those parameters there are no subjects "+\
            "available in this study, make the selectors files are correct"

    return chosen_samples, final_biom_table
