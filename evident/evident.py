__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011-2012, Evident"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "0.1.dev"
__maintainer__ = ["Antonio Gonzalez Pena"]
__email__ = "antgonza@gmail.com"
__status__ = "Development"

from numpy import Inf, array, abs

"""
HSV colors taken from qiime.colors
Reformated data:
"""
data_color_hsv = [
#(0.0,0.0,0.20),
(0.0,0.89,0.89),
(0.240,0.89,0.89),
#(0.28,0.98,0.95),
(0.120,0.89,0.50),
(0.302,0.73,0.57),
(0.60,0.89,0.89),
(0.184,0.49,0.96),
(0.333,0.37,0.96),
(0.178,0.42,0.63),
(0.36,0.89,0.42),
(0.0,0.0,0.50),
(0.123,0.99,0.96),
(0.14,0.51,0.97),
(0.211,0.42,0.85),
(0.32,0.46,0.99),
(0.142,0.36,0.79),
(0.269,0.29,0.75),
(0.56,0.40,0.89),
(0.303,0.89,0.24),
(0.0,0.0,0.75),
#(0.192,0.89,0.24),
(0.325,0.89,0.93),
(0.197,0.89,0.89),
(0.271,0.43,0.36),
(0.33,0.45,0.77),
(0.60,0.89,0.50),
(0.264,0.75,0.89),
#(0.60,0.66,0.75),
#(0.213,0.45,0.77),
(0.348,0.31,0.74),
(0.180,0.89,0.50),
#(0.60,0.89,0.28),
(0.0,0.89,0.50),
(0.81,0.89,0.26),
#(0.240,0.89,0.41),
(0.26,0.89,0.65),
#(0.25,0.89,0.20),
(0.17,0.89,0.63),
(0.272,0.89,0.44)
]

def make_pcoa_plot(ellipses, mapping, pcoalabels):
    """ creates the webGL output to generate the PCoA plots
    
    Input:
    samples { 
        'sample_1': { 
               'center': array of n dimensions, original pcoa axis, 
               'radii': array of n dicts {
                           'max': max value of randomizations
                           'avg': avg value of randomizations
                           'min': min value of randomizations
               }
        }
    }
    
    Output:
        webGL text
    """
    
    result = "var headers = ["
    for h in mapping[1]:
        result += "\""+h+"\","
    result = result[:-1]+"]\n"
    
    result += "var mapping = {"
    for l in mapping[0]:
        result += "\""+l[0]+"\": ["
        for i in range(len(l)):
            result += "\""+l[i]+"\","
        result = result[:-1]+"],"
    result = result[:-1]+"}\n"
    
    # Generating points
    min_coords = array([ Inf, Inf, Inf])
    max_coords = array([Inf*-1, Inf*-1, Inf*-1])
    
    result += 'var ellipses = new Array()\n'
    result += 'var points = new Array()\n'
    for i, (sample,values) in enumerate(ellipses.items()):
        result += "ellipses['%s'] = { 'name': '%s', 'color': %d, 'width': %f, 'height': %f, 'length': %f , 'x': %f, 'y': %f, 'z': %f }\n" % \
            (sample, sample, 0, values['axes_radii'][0], values['axes_radii'][1], values['axes_radii'][2], values['center'][0], values['center'][1], values['center'][2])
        result += "points['%s'] = { 'name': '%s', 'color': %d, 'x': %f, 'y': %f, 'z': %f }\n" % \
            (sample, sample, 0, values['center'][0], values['center'][1], values['center'][2])
        
        if min_coords[0]>values['center'][0]-values['axes_radii'][0]: min_coords[0]=values['center'][0]-values['axes_radii'][0]
        if min_coords[1]>values['center'][1]-values['axes_radii'][1]: min_coords[1]=values['center'][1]-values['axes_radii'][1]
        if min_coords[2]>values['center'][2]-values['axes_radii'][2]: min_coords[2]=values['center'][2]-values['axes_radii'][2]
        
        if max_coords[0]<values['center'][0]+values['axes_radii'][0]: max_coords[0]=values['center'][0]+values['axes_radii'][0]
        if max_coords[1]<values['center'][1]+values['axes_radii'][1]: max_coords[1]=values['center'][1]+values['axes_radii'][1]
        if max_coords[2]<values['center'][2]+values['axes_radii'][2]: max_coords[2]=values['center'][2]+values['axes_radii'][2]
    
    result += 'var segments = 16, rings = 16, radius = %f\n' % ((min_coords[0]-max_coords[0])*.02)
    result += 'var xaxislength = %f\n'%(max_coords[0]+abs(min_coords[0]))
    result += 'var yaxislength = %f\n'%(max_coords[1]+abs(min_coords[1]))
    result += 'var zaxislength = %f\n'%(max_coords[2]+abs(min_coords[2]))
    result += 'var min_x = %f\n' % min_coords[0]
    result += 'var min_y = %f\n' % min_coords[1]
    result += 'var min_z = %f\n' % min_coords[2]
    result += 'var max_x = %f\n' % max_coords[0]
    result += 'var max_y = %f\n' % max_coords[1]
    result += 'var max_z = %f\n' % max_coords[2]
    
    result += "max = %f\n" % (max([abs(max_coords).max(), abs(min_coords).max()])*1.15)
    # result += "var s = \"%s\"\n" %(pcoalabels)
    result += "pc1 = %.0f\n" %(pcoalabels[0])
    result += "pc2 = %.0f\n" %(pcoalabels[1])
    result += "pc3 = %.0f\n" %(pcoalabels[2])
    return result
