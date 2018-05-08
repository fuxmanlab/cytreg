# -*- coding: utf-8 -*-

def convert_dict_to_utf8(the_dict):
    ''' Convert all keys and values (stored as lists)
        to UTF-8 encoding '''
    output = {}
    for key in the_dict:
        output[key.decode('utf-8').encode('utf-8')] = [x.decode('utf-8').encode('utf-8') for x in the_dict[key]]
    return output


def read_tf_list(filename):
    '''  Read the TF list file into a dictionary of tFName:[tfAlias].
         Convert all strings to UTF-8 encoding '''
    tf_list = []
    with open(filename) as f:
        for line in f.readlines():
            tf_list.append(tuple( (val.strip() for val in line.split(' ',1) ) ))
    # Initialize the dictionary from all the unique tfNames
    tfNames = set([x[0] for x in tf_list])
    tfDict = {key:[] for key in tfNames}
    # Loop thru tf_list and append all tfAlias values to the appropriate key
    for tf in tf_list:
        tfDict[tf[0]].append(tf[1])
    return tfDict



# definitions for NIH EUtils scripts

# cytokines
# NOTE:!!!  This must include beta and alpha characaters where appropriate for proper abstract matching!
# example:  'TGF-β'
# Force cyts into a complete 8-bit UTF-8 encoding
#cyts = convert_dict_to_utf8(cyts)

# assays
# In order of most to least frequently found - speeds up matching in post_process_pmids.py
assays = ['Chromatin+Immunoprecipitation','Chromatin+Immuno-precipitation','Chromatin-Immunoprecipitation','Chromatin-Immuno-precipitation','Chromatinimmunoprecipitation','Ch-IP','EMSA','electrophoretic+mobility+shift+assay','mobility+shift+electrophoresis','gel+shift+assay','gel+mobility+shift+assay','band+shift+assay','gel+retardation+assay','Reporter+Assay','Functional+Assay','luciferase+assay','promoter+assay','promoter+activity+assay','Transformation+Assay','Transfection+assay','gene+expression+assay','GUS+Assay','Functional-Assay','luciferase-assay','promoter-assay','promoter-activity+assay','Transformation-Assay','Transfection-assay','gene-expression+assay','GUS-Assay','Reporter-Assay','promoter-activity-assay','gene-expression-assay']
# Convert assays as well, just to be sure.
assays = [a.encode('utf-8') for a in assays]

#helper_body_terms = ['enhancer','promoter','upstream']
#helper_body_terms = [a.encode('utf-8') for a in helper_body_terms]

#helper_sentence = ['induc','modulat','promot','enhance','regulat','bind','activat','repress','mediat']
#helper_sentence = [a.encode('utf-8') for a in helper_sentence]
