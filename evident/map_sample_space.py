#!/usr/bin/env python
# File created on 07 Mar 2012
#from __future__ import division

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2012, Evident"
__credits__ = ["Will Van Treuren", "Antonio Gonzalez Pena",
"Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"


from numpy import searchsorted
from qiime.filter import filter_mapping_file
from qiime.util import compute_seqs_per_library_stats

def get_sorted_counts_per_sample(biom_table, reverse=False):
    """gets a sorted list of sequences per sample from min to max
    
    inputs:
    biom_table: biom table object
    revers: reverse the ordering value i. e. from max to min
    
    outputs:
    sorted_counts_per_sample: list of tuples sorted on first element which
    gives [(seqs/sample, sampleId)... ]
    """

    sample_counts = compute_seqs_per_library_stats(biom_table)[4]
    
    sorted_counts_per_sample = [(v,k) for k,v in sample_counts.items()]
    sorted_counts_per_sample.sort()
    
    if reverse:
        sorted_counts_per_sample.reverse()

    return sorted_counts_per_sample

def make_selectors(counts_per_sample, minimum, mapping_file_tuple,
                    subject_header_name, verbose=False):
    """make the four column string needed to print in the selectors file

    Inputs:
    counts_per_sample: a sorted list of tuples with the sample identifier and
    the number of sequences.
    minimum: minimum number of sequences considered to be a valid state.
    mapping_file_tuple: a tuple with the data of a mapping file and the headers.
    subject_header_name: string identifying the name of the column in the 
    mapping file that represents a unique subject.

    Output:
    result: four columns string corresponding to number of sequences, subjects,
    number of samples and metadata fields.
    """

    # unwrap the mapping file
    mapping_data = mapping_file_tuple[0]
    mapping_headers = mapping_file_tuple[1]

    seqs_per_sample = [t[0] for t in counts_per_sample]

    head_val = None
    subj_val = None
    samp_sub = None
    results = []

    depth = -1
    samples_per_subject = {}

    # store the index for convenience
    subject_index = mapping_headers.index(subject_header_name)
    list_of_subjects = [line[subject_index] for line in mapping_data]

    # initialize the samples_per_subject dictionary with as many keys as
    # subjects and values equal to the minimum number of samples among them
    for unique_subject in list(set(list_of_subjects)):
        samples_per_subject[unique_subject] = list_of_subjects.count(unique_subject)
    least_number_of_samples = min(samples_per_subject.values())
    for key, value in samples_per_subject.iteritems():
        samples_per_subject[key] = least_number_of_samples

    for sequences_per_sample_tuple in counts_per_sample:

        # there's no need to iterate if the minimum rarefaction depth is not met
        # or if the depth is the same as the previous depth, this would mean a 
        # repeated row in the output line with the same values
        if  sequences_per_sample_tuple[0] < minimum or \
            sequences_per_sample_tuple[0] == depth:
            continue

        if verbose:
            print 'Samples per subject: {0} @ depth: {1}'\
                .format(samples_per_subject, depth)

        # Some samples are not in the mapping file just print those out
        sample_id = sequences_per_sample_tuple[1]
        try:
            current_subject = [line[subject_index] for line in mapping_data if line[0] == sample_id][0]
        except IndexError:
            print 'Sample Id: {0} is not in the mapping file'.format(sample_id)
            continue

        # extract convenience data for ease of use
        depth = sequences_per_sample_tuple[0]
        remaining_ids = [_tuple[1] for _tuple in counts_per_sample if _tuple[0] >= depth]

        filtered_headers, filtered_data = filter_mapping_file(mapping_data,\
            mapping_headers, remaining_ids, include_repeat_cols=False)

        # Breaking when there are no subjects/individuals left
        if subject_header_name not in filtered_headers:
            break

        # numbers to be written in the selectors file
        number_of_subjects = len(samples_per_subject.keys())
        number_of_samples = min(samples_per_subject.values())

        if number_of_subjects*number_of_samples < 3:
            continue

        # format the output 
        if not subj_val and not head_val and not samp_sub:
            results.append('%d\t%d\t%d\t%s' % (int(depth), number_of_subjects,\
                number_of_samples, ','.join(filtered_headers[1:-1])))
            subj_val = number_of_subjects
            head_val = filtered_headers
            samp_sub = number_of_samples
        else:
            if head_val!=filtered_headers:
                results.append('%d\t%d\t%d\t%s'%(int(depth),number_of_subjects,\
                    number_of_samples, ','.join(filtered_headers[1:-1])))
                head_val = filtered_headers
            elif samp_sub!=number_of_samples:
                results.append('%d\t%d\t%d\tNone'% (int(depth),\
                    number_of_subjects, number_of_samples))
                samp_sub = number_of_samples
            elif subj_val!=number_of_subjects:
                results.append('%d\t%d\t%d\tNone' % (int(depth),\
                    number_of_subjects, number_of_samples))
                subj_val = number_of_subjects

        # remove the current processed sample and if needed, remove the subject
        try:
            samples_per_subject[current_subject] -= 1
            if samples_per_subject[current_subject] == 0:
                del samples_per_subject[current_subject]
        except:
            pass

    return results
