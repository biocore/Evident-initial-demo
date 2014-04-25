__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011-2012, Evident"
__credits__ = ["Antonio Gonzalez Pena", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "0.1.dev"
__maintainer__ = ["Antonio Gonzalez Pena"]
__email__ = "antgonza@gmail.com"
__status__ = "Development"


from flask import Flask, render_template, session, request
app = Flask(__name__, template_folder='www', static_folder='www')
app.secret_key = 'evident_secret_key'

from evident.parse import load_studies

from qiime.filter import filter_mapping_file
from qiime.format import format_biom_table
from qiime.parse import parse_mapping_file, parse_newick

from biom.parse import parse_biom_table
from biom.exception import TableException

# from evident.error import raiseApacheError
from evident.subsampling import select_samples
from evident.rarefaction import generate_alpha_rarefaction_data_from_point_in_omega

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


if __name__ == '__main__':
    app.run(debug=True)
