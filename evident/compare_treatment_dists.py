#!/usr/bin/env python
# File created on 29 Jan 2013
from __future__ import division

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2013, The Evident Project"
__credits__ = ["Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"

from numpy import arange, eye, nan, zeros
from qiime.beta_diversity import single_object_beta
from qiime.parse import parse_distmat
from numpy.ma import masked_array
from numpy.ma.core import MaskedConstant
from itertools import combinations

def between_treatments_dist(group_t1, group_t2, marginal_ids, distmat):
    """Calculate avg dist, se between treatment 1 and treatment 2 sample groups.
    Notes:
     group_t1 - list of ids found in marginal ids.
     group_t2 - list of ids found in marginal ids. 
     marginal_ids - list of ids, identifies the cols/rows of the distmat.
     distmat - symmetric hollow array, dist between samples.
     This function calculates the average distance between all samples in 
     group_t1 and all samples in group_t2. 
    """
    if type(marginal_ids)!=list:
        marginal_ids = list(marginal_ids)
    g1 = [marginal_ids.index(i) for i in group_t1]
    g2 = [marginal_ids.index(i) for i in group_t2]
    d = distmat.take(g1,0).take(g2,1)
    return d.mean(), d.std()/d.size

def treatment_dist(group, marginal_ids, distmat):
    """Calc avg dist, se between group and all other samples in distmat.
    Notes:
     group - list, group of samples to calc intersample dist of.
     distmat - symmetric hollow 2d numpy array, dist between samples listed in 
     marginal_ids in order of marginal_ids.
     marginal_ids - list, identifies the cols/rows of the distmat.
     formula:
     B = set(marginal_ids)-set(group) # ids in the distmat that arent in group
     s = sample id in group
     k = sample id in B
     summation(dist(s,k))/[len(B)*len(group)]
     excludes dist to self D(i,i) from averages.
    """
    if type(marginal_ids)!=list:
        marginal_ids = list(marginal_ids)
    mdi = list(arange(len(marginal_ids))) # marginal indices indexes in dm
    group_indices = [marginal_ids.index(i) for i in group]
    r = distmat[group_indices,:].take(list(set(mdi)-set(group_indices)),1)
    return r.mean(), r.std()/r.size

def within_treatment_dist(group, marginal_ids, distmat):
    """Calc avg dist in group.
    Notes:
     group - list, group of samples to calc intersample dist of.
     distmat - symmetric hollow 2d numpy array, dist between samples listed in 
     marginal_ids in order of marginal_ids.
     marginal_ids - list, identifies the cols/rows of the distmat.
     excludes dist to self D(i,i) from averages.
    """
    if type(marginal_ids)!=list:
        marginal_ids = list(marginal_ids)
    group_indices = [marginal_ids.index(i) for i in group]
    r = distmat.take(group_indices,0).take(group_indices,1)
    num_comps = float(r.size - r.shape[0]) # r is square, r.shape[0]=r.shape[1]
    # unsure about the se calculation -- we are trying to calculate the
    # standard error of this hollow symmetric distmat. however, because of its
    # symmetric nature, i am not sure if i should divide by num_comps or by
    # just r.shape[0](r.shape[0]-1)/2 since i am counting redundant information
    # for reduction of the se in the first case.
    se = masked_array(r, eye(r.shape[0])).std()/num_comps
    if type(se) == MaskedConstant: # array was size one, no se
        se = nan 
    return r.sum()/num_comps, se

def treatment_covering(sids, category, mf):
    """Return dict of treatment to sample id list for sample id in sids.
    Notes:
     sids - list coming from session['chosen_samples'], sample ids.
     e.g [sam1,sam4,sam7,sam3,sam6,...]
     category - str, category under the mf that samples are to be comapred on.
     mf - parsed mapping file, dict of sample_id:metadata.
     output - {treatment:[list of samples in sids who have that value]}
    """
    # Yoshiki found a more elegant way than chain:
    # list(set([mf[i][category] for i in list(chain.from_iterable(sids))]))
    # the new version of evident appears not to five sids as a single list
    # flat_sids = sum(sids,[])
    # treatments = list(set([mf[i][category] for i in flat_sids]))
    treatments = list(set([mf[i][category] for i in sids]))
    return {t:[i for i in sids if mf[i][category]==t] for t in treatments}

def compare_treatment_dists(chosen_samples, category, mf, bt, m, tr):
    """Calculate avg between, within, and to-all distances for chosen_samples.
    Notes: 
     chosen_samples is a list of lists of ids that collectively have some amount
     of different values under category in the mapping file. these samples will
     br grouped by the value they have and then these groupings will be 
     compared. the between distance is the all the pairwise distances between 
     the groupings. the within distance is the distance between the samples in a 
     single group. the to-all distance is the distance from the group to all
     other samples in the distmat.
    Inputs:
     chosen_samples - list of ids. e.g. [sam1,sam7,sam3,sam6,..]
     category - str, field in mf.
     mf - parsed mapping file, dict of sample_id:metadata.
     bt - biom table containing at least all samples contained in the mf.
     m - str, metric to used for beta diversity calculation. 
     tr - tree object, containing at least all nodes in bt.
    Output:
     A list of marginals that are the treatments of the groups, i.e. ['HF','LF']
     bt_wi_m - a 2d upper triangular array that has the average distances
     between treatment groups (or in the case of the main diagonal, the average
     within treatment/group distance).
     bt_wi_se - the standard errors for bt_wi_m.
     ta_m_se - 2d array with number of treatments/groups rows, and 2 cols where 
     the first col is the average distance between that treatment and all others
     and the second col is the se. 
    """
    dm = single_object_beta(bt, m, tr) #make the sample-sample distance matrix
    samples, data = parse_distmat(dm) #parse dm which is list of strs
    tc = treatment_covering(chosen_samples, category, mf)
    output_marginals = tc.keys()
    # make 3 arrays for output, between-within means, between-within ses, 
    # to-all means and ses,
    bt_wi_m = zeros((len(output_marginals),len(output_marginals)))
    bt_wi_se = zeros((len(output_marginals),len(output_marginals)))
    ta_m_se = zeros((len(output_marginals),2))
    for i,t in enumerate(output_marginals): # calculate within and to-all
        ta_m_se[i][0], ta_m_se[i][1] = treatment_dist(tc[t], samples, data)
        bt_wi_m[i][i], bt_wi_se[i][i] = within_treatment_dist(tc[t], samples,
            data)
    for t1, t2 in combinations(output_marginals, 2): #calculate between dists
        t1_ind = output_marginals.index(t1)
        t2_ind = output_marginals.index(t2)
        bt_wi_m[t1_ind][t2_ind], bt_wi_se[t1_ind][t2_ind] = \
            between_treatments_dist(tc[t1], tc[t2], samples, data)
    return output_marginals, bt_wi_m, bt_wi_se, ta_m_se




