#!/usr/bin/env python

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011-2012, Evident"
__credits__ = ["Antonio Gonzalez Pena", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "0.1.dev"
__maintainer__ = ["Antonio Gonzalez Pena"]
__email__ = "antgonza@gmail.com"
__status__ = "Development"

# prevents runtime error with matplotlib, hence this statement must preceed
# any interaction with that module to avoid unwanted results
import os, tempfile, logging
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

from StringIO import StringIO

from flask import (Flask, render_template, request, make_response, session,
                   url_for)
from qiime.filter import filter_mapping_file
from qiime.format import format_biom_table
from qiime.parse import parse_mapping_file, parse_newick, parse_coords
from qiime.make_rarefaction_plots import make_averages
from biom.parse import parse_biom_table
from biom.exception import TableException

from evident.demo import demo
from evident.subsampling import select_samples
from evident.rarefaction import (build_color_preferences,
    generate_alpha_rarefaction_data_from_point_in_omega)
from evident.pcoa import (generate_pcoa_cloud_from_point_in_omega,
                          make_pcoa_plot)
from evident.parse import load_studies
from evident.session import SqliteSessionInterface


# set up the Evident application
app = Flask(__name__, template_folder='www', static_folder='www')


@app.route('/')
def init():
    return render_template('index.html')


@app.route('/_load_study_data')
def load_study_data():
    return load_studies()


@app.route('/_subsample')
def subsample():
    # loading the gg tree can take a while; load it when alpha is selected
    session['tree_object'] = None

    session['study'] = request.args.get('study', None, type=str)
    session['sequences'] = request.args.get('sequences', None, type=int)
    session['demo'] = request.args.get('demo', None, type=str)
    session['iterations'] = request.args.get('iterations', None, type=int)
    subjects = request.args.get('subjects', None, type=int)
    samples = request.args.get('samples', None, type=int)

    session['pcoa'] = request.args.get('pcoa', None, type=str)
    session['alpha_stddev'] = request.args.get('alpha_stddev', None, type=str)
    session['alpha_stderr'] = request.args.get('alpha_stderr', None, type=str)

    # Validating arguments
    study = {}
    for i,line in enumerate(open('data/studies.txt')):
        line = line.strip().split('\t')
        if len(line)<1 or i==0 or line[0].startswith('#'):
            continue
        if str(session['study']) == str(line[1]):
            study['name'] = line[0]
            session['filename'] = line[1]
            study['subject_column'] = line[2]

    # If the study couldn't be found or parsed
    if study == {}:
        raise ValueError('<b>Study <u>%s</u> does not exist.</b>' % form_input['study'])

    # Creating full paths to files
    mapping_fp = 'data/' + session['filename'] + '_map.txt'
    biom_fp = 'data/' + session['filename'] + '.biom'
    alpha_fp = 'data/' + session['filename'] + '_alpha.html'

    if session['demo']!="true":
        map_data, headers, comments = parse_mapping_file(open(mapping_fp, 'U'))
        biom_table = parse_biom_table(open(biom_fp, 'U'))

        try:
            # select only samples that meet ther criteria
            chosen_samples, filtered_biom_table = select_samples(map_data, headers, biom_table, session['sequences'], study['subject_column'], subjects, samples)

            # check if we have enough samples to display a PCoA plot
            if len(chosen_samples) < 3:
                raise ValueError('<b>At least <u>three</u> data-points are needed: try changing the values of Subjects or Samples per Subject.</b>')

            session['chosen_samples'] = chosen_samples
            session['filtered_biom_table'] = format_biom_table(filtered_biom_table)
        except TableException:
            raise ValueError('<b>There are <u>not enough</u> subjects in this study: check your BIOM file and Mapping File have equal number of samples.</b>')

        # get the non redundant columns of the mapping file
        interesting_headers, interesting_data = filter_mapping_file(map_data, headers, chosen_samples)

        # specific to alpha diversity visualizations
        if session['alpha_stddev'] or session['alpha_stderr']:
            tree_fp = 'data/gg_97_otus_4feb2011.tre'
            session['tree_object'] = parse_newick(open(tree_fp))

            session['alpha_rarefaction_data'] = generate_alpha_rarefaction_data_from_point_in_omega(biom_object=filtered_biom_table, metrics=['observed_species','Chao1','PD_whole_tree'], sequences=session['sequences'], iterations=session['iterations'], tree_object=session['tree_object'])

    # demo
    else:
        interesting_data, interesting_headers, comments = parse_mapping_file(open(mapping_fp,'U'))
        session['interesting_headers'] = interesting_headers
        session['interesting_data'] = interesting_data

    session['mapping_file_tuple'] = (interesting_data, interesting_headers)

    return session['demo']

@app.route('/_results')
def results():
    max_iterations = 10
    min_sequences  = 10

    # Reading / parsing some values
    iterations = session['iterations']
    print session['pcoa'], session['alpha_stderr'], session['alpha_stddev']
    viz = request.args.get('viz', None, type=str)

    # Validating that number of iterations
    if max_iterations<iterations:
        raise ValueError('<b>Sorry we can not process more than %s iterations.</b>' % (max_iterations))

    # Creating full paths to files
    mapping_fp = 'data/' + session['filename'] + '_map.txt'
    pcoa_fp = 'data/' + session['filename'] + '_unweighted_unifrac_pc.txt'
    alpha_stddev_fp = 'data/' + session['filename'] + '_alpha_stddev.html'
    alpha_stderr_fp = 'data/' + session['filename'] + '_alpha_stderr.html'

    # not demo, live user interaction
    if session['demo']!="true":

        # the gg tree is only loaded if it wasn't loaded in lib.psp
        if session['tree_object'] == None:
            tree_fp = 'data/gg_97_otus_4feb2011.tre'
            tree_object = parse_newick(open(tree_fp))
        else:
            tree_object = session['tree_object']

        mapping_file_tuple = session['mapping_file_tuple']
        filtered_biom_table = parse_biom_table(StringIO(session['filtered_biom_table']))

        # principal coordinates analysis plots
        if viz=='pcoa':
            # try:
            webgl_string = generate_pcoa_cloud_from_point_in_omega(
                    map_headers=mapping_file_tuple[1],
                    map_data=mapping_file_tuple[0],
                    biom_object=filtered_biom_table, metric='unifrac',
                    sequences=session['sequences'], iterations=iterations, axes=3,
                    tree_object=tree_object)
            location =  url_for('static', filename='emperor_required_resources')
            return (webgl_string).replace('emperor_required_resources', location)

        # alpha rarefaction plots
        elif viz=='alpha_stddev' or viz=='alpha_stderr':
            min_depth = 10
            if session['sequences']==min_depth:
                raise ValueError('The min number of sequences for alpha diversity is 15, so you need to select at least this number')

            # create all the coloring data for the alpha rarefaction plots
            prefs, data, background_color, label_color = build_color_preferences(
                mapping_file_tuple)

            # make the rarefaction plots and rander the results on screen
            html_string = make_averages(prefs, data, background_color, label_color,
                session['alpha_rarefaction_data'], 'dummy_fp', 75, 'png', None,
                False, viz[6:], 'memory')
            return (html_string)

        else:
            raise ValueError('Visualization <u>%s</u> does not exist' % viz)

    # demos
    else:
        # principal coordinates analysis plots
        if viz=='pcoa':
            map_data, map_headers, comments = parse_mapping_file(open(mapping_fp,'U'))
            pcoa_headers, pcoa_values, eigenvalues, coords_pct = parse_coords(open(pcoa_fp,'U'))

            webgl_string = make_pcoa_plot(pcoa_headers, pcoa_values, eigenvalues, coords_pct,
                map_headers, map_data)

            # to link to local resources and avoid forcing internet access
            location =  url_for('static', filename='emperor_required_resources')
            return webgl_string.replace('emperor_required_resources', location)

        # alpha rarefaction plots
        elif viz=='alpha_stddev' or viz=='alpha_stderr':
            if viz == 'alpha_stddev':
                html_string = '\n'.join(open(alpha_stddev_fp,'U').readlines())
            elif viz == 'alpha_stderr':
                html_string = '\n'.join(open(alpha_stderr_fp,'U').readlines())

            return (html_string)
        else:
            raise ValueError('Visualization <u>%s</u> does not exist' % viz)


if __name__ == '__main__':
    # Starting sessions
    # You can change the session storage to wherever you want
    app.session_interface = SqliteSessionInterface('/tmp/evident/')
    app.run(debug=True, port=8888)
