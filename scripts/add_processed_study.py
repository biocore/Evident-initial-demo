#!/usr/bin/env python
# File created on 01 Apr 2013
# This script is based on: beta_diversity_through_plots.py, alpha_rarefaction.py
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2012, Evident"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = ".9-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"


from qiime.util import parse_command_line_parameters, make_option, get_options_lookup
from os.path import exists, join
from os import listdir
from re import sub
from shutil import copyfile, move
from datetime import datetime
options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Add a study to Evident"
script_info['script_description'] ="""This script takes a folder created by \
process_new_study.py and adds the study into Evident. It validates the presence of the\
files but not their validity.

Note: If you use the name of an existing study, the script will replace this.
"""

script_info['script_usage'] = [("Add a study using the defaults","",
'%prog -i crawford_preprocessed/ -s "Crawford et al. 2009 - mice fasting"')]
script_info['script_usage'].append(("Add a study and forcing the output","",
'%prog -i crawford_preprocessed/ -s "Crawford et al. 2009 - mice fasting" -f'))

script_info['output_description']= """The script doesn't have a direct output but will \
add a new study to Evident"""
script_info['required_options'] = [\
 make_option('-i','--input_path',type='existing_path',
            help='the input folder, the result of process_new_study.py [REQUIRED]'),
 make_option('-s','--study_name',type='string', help='the name of the study [REQUIRED]'),           
]
script_info['optional_options'] = [\
 make_option('-o','--output_path',type='existing_path',
            help='the output folder, the data folder for Evident [default: %default]',
            default='/evident/data/'),
 make_option('-f','--force_overwrite',action='store_true',
            help='force the overwrite of output files [default: %default]',
            default=False),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    # converting the name of the study to a valid file name
    study_prefix = sub(r'[^\w]+','_',opts.study_name).lower()
    study_name = join(opts.output_path, study_prefix)
    
    # validating the preexistance of the output files
    if (exists(study_name + '.biom') or exists(study_name + '_alpha_stderr.html') or \
       exists(study_name + '_selectors.txt') or \
       exists(study_name + '_alpha_stddev.html') or exists(study_name + '_map.txt') or \
       exists(study_name + '_unweighted_unifrac_pc.txt')) and not opts.force_overwrite:
       raise IOError, 'The output file(s) exist, either change the name of the study ' +\
          'or use -f'
    
    # validating existance of input files
    biom_file = join(opts.input_path, 'raw.biom')
    if not exists(biom_file):
        raise IOError, "Couldn't find the raw biom file in the input folder"
        
    alpha_stderr = join(opts.input_path, 'alpha/alpha_rarefaction_plots_stderr/rarefaction_plots.html')
    if not exists(alpha_stderr):
        raise IOError, "Couldn't find an alpha stderr file in the input folder"
    
    alpha_stddev = join(opts.input_path, 'alpha/alpha_rarefaction_plots_stddev/rarefaction_plots.html')
    if not exists(alpha_stddev):
        raise IOError, "Couldn't find an alpha stddev file in the input folder"
    
    selectors = join(opts.input_path, 'selectors.txt')
    if not exists(selectors):
        raise IOError, "Couldn't find a selectors file in the input folder"
    
    mapping_file = join(opts.input_path, 'mapping_file.txt')
    if not exists(mapping_file):
        raise IOError, "Couldn't find a mapping file in the input folder"
    
    unweighted_unifrac_pc = join(opts.input_path, 'unweighted_unifrac_pc.txt')
    if not exists(unweighted_unifrac_pc):
        raise IOError, "Couldn't find a unweighted unifrac pc file in the input folder"
    
    
    # coping files to destination 
    copyfile(biom_file, join(opts.output_path, study_name + '.biom'))
    copyfile(alpha_stderr, join(opts.output_path, study_name + '_alpha_stderr.html'))
    copyfile(alpha_stddev, join(opts.output_path, study_name + '_alpha_stddev.html'))
    copyfile(selectors, join(opts.output_path, study_name + '_selectors.txt'))
    copyfile(mapping_file, join(opts.output_path, study_name + '_map.txt'))
    copyfile(unweighted_unifrac_pc, join(opts.output_path, study_name + '_unweighted_unifrac_pc.txt'))
    
    # modifying the studies.txt file
    lines = open(join(opts.input_path, 'study_preferences.txt'),'U').read().split('\n')
    rarefied = lines[0]
    ind_id = lines[1]
    
    content = {}
    for i,line in enumerate(open(join(opts.output_path, 'studies.txt'),'U')):
        line = line.strip().split('\t')
        if i==0:
            headers = line
        else:
            content[line[1]] = line
    
    new_file = 'studies_%s.txt' % datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    new_study = [opts.study_name, study_prefix, ind_id, rarefied]
    content[study_prefix] = new_study
    
    results = ['\t'.join(headers)]
    for k,v in content.items():
        results.append('\t'.join(v))
        
    move(join(opts.output_path, 'studies.txt'),join(opts.output_path, new_file))
    open(join(opts.output_path, 'studies.txt'),'w').write('\n'.join(results))
    
    
if __name__ == "__main__":
    main()
