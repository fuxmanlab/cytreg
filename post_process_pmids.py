# -*- coding: utf-8 -*-

#  Sample run command:
# python post_process_pmids_0.4.py -xml xmlDir -tf_list ../tf_list.txt -match xml -ofile test.out -a

# Recommended: use Intel optimized python for 20x faster BeautifulSoup!!!


############ Version 0.4  4/10/2017
#
#  Added argparse arguments
#  Added substring matching option for assays
#
############ Version 0.3    4/5/2017
#   Changed string matching to regex to properly match words.
#    
#
############ Version 0.1     3/22/207
#   Initial version
#
############

# Import local data definitions
from nih_eutils_defs import *


# Process the data in parallel...
import concurrent.futures

import re  # regular expression library3
import argparse  # command line arguments
from bs4 import BeautifulSoup
import os
import datetime
import fnmatch
import io
import codecs
import traceback
import sys
UTF8Writer = codecs.getwriter('utf8')
sys.stdout = UTF8Writer(sys.stdout)

 
class PmidRecord(object):
    def __init__(self,xml_soup, ignorecase=True):
        ''' xml_soup is a BeautifulSoup article.
            Load all the text as 8-bit strings in UTF-8 encoding. '''
        # The default behaviour of BSoup is to encode in UTF-8
        self.xml = xml_soup.get_text().encode('utf-8')
        self.pmid = int(xml_soup.find('pmid').get_text())
        self.abstract = xml_soup.find('abstract')
        try:
            if self.abstract: # Sometimes these are missing!
                # No parts just use the whole abstract
                if xml_soup.find('abstract').findChildren('abstracttext'):
                    abstract_parts = xml_soup.find('abstract').findChildren('abstracttext')
                    abstract = io.StringIO(u'')
                    for part in abstract_parts:
                        tmp = abstract.write('%s\n\n' % part.get_text())
                    self.abstract = abstract.getvalue().encode('utf-8')
                else:
                    self.abstract = xml_soup.find('abstract').getvalue().encode('utf-8')
        except: # Just give up.
            self.abstract = None
        self.title = xml_soup.find('articletitle').get_text().encode('utf-8')
        # Store case matching
        self.ignorecase = ignorecase
        self._flags = re.UNICODE | re.IGNORECASE if self.ignorecase else re.UNICODE


    def set_match(self,selection):
        ''' Pick a match method based on a keyword'''
        if selection == 'xml':
            return self.match_xml
        elif selection == 'abstract':
            return self.match_abstract
        elif selection == 'title':
            return self.match_title
        raise Exception('Match options are: xml, abstract, title')

    # The "x in y" syntax is the fastest way to do string matching in Python
    def match_xml(self,term, strict=True):
        ''' Match one term in the XML field '''
        return self._match_term_to_field(term,self.xml, strict)
        
    def match_abstract(self,term, strict=True):
        ''' Match one term in the abstract field '''
        if self.abstract:   
          return self._match_term_to_field(term,self.abstract, strict)
        return False # doesn't exist!

    def match_title(self,term, strict=True):
        '''Match one term in the abstract field '''
        return self._match_term_to_field(term,self.title, strict)

    def _match_term_to_field(self,term,field, strict=True):
        ''' Match a search term to a field in the object '''
        if strict:
            regex = re.compile(r'\b%s\b' % re.escape(term),flags = self._flags)
            return regex.search(field) is not None  
        ''' Otherwise just do a substring match.  The Python 'in' syntax
            is faster than another regex '''
        return term.lower() in field.lower() >= 0


class ProcessXmlDirectory(object):
    ''' This class processes a directory of PubMed XML downloads and 
        writes matches to a file. '''
    def __init__(self,directory,output_file,match_type,assays,cyt_dict,tf_dict, strict_assay=True):
        self.directory = directory
        self.output_file = output_file
        self.match_type = match_type
        self.assays = assays
        self.cyt_dict = cyt_dict
        self.tf_dict = tf_dict
        self.strict_assay = strict_assay
        
    def process_xml_file(self,xml_file):
        ''' Process an XML file and return matches for assays, cytokines,
            and TFs.   This is written as a generator for an external
            for loop to handle.
            
            Inputs:  xml_file - XML filename
                     match_type - what to match to: xml,abstract,title
                     assays - list of assays
                     cyts_dict - cytokine dictionary
                     tf_dict - TF dictionary 
            Output:  String in the format     cyt_name,cyt_alias,tf_name,tf_alias,assay,rec.pmid
        '''
        # Turn the XML file to soup.  This is slow, ~20 seconds per file with regular Python.
        # It is much faster, ~1.5 seconds, with the Intel optimized Python.
        try:
            with open(xml_file) as f:
                soup = BeautifulSoup(f,'lxml')
            # Locate all the articles
            articles = soup.find('pubmedarticleset').findChildren('pubmedarticle')
            # Loop over all articles
            if articles:
                for xml_soup in articles:
                    rec_s = datetime.datetime.now()
                    rec = PmidRecord(xml_soup, ignorecase=True)
                    match = rec.set_match(self.match_type)
                    for assay in self.assays:
                        if match(assay, strict=self.strict_assay):
                            for cyt_name in self.cyt_dict:
                                for cyt_alias in cyts[cyt_name]:
                                    if match(cyt_alias):
                                        for tf_name in self.tf_dict:
                                            for tf_alias in tf_dict[tf_name]:
                                                if match(tf_alias):
                                                    yield '%s, %s, %s, %s, %s, %s' % (cyt_name,cyt_alias,tf_name,tf_alias,assay,rec.pmid)                 
        except Exception, e:
            # Something went wrong.  Print out the filename with its path and the exception.
            # Don't yield anything so that the output file doesn't get misformatted.
            msg = 'Error processing file: %s\n%s' % (xml_file,e)
            sys.stderr.write('%s\n' % msg)
            traceback.print_exc(file=sys.stderr)
            sys.stderr.flush()
                        
    def process_xml_directory(self):
        ''' Takes a directory containing the xml output files from
            PubMed.  Loops over them and returns the output of 
            process_xml_file.  Runs as a generator. To use this a double
            for loop is needed - one to loop over the files, the other to
            loop over the output per file.
        '''
        for filename in os.listdir(self.directory):
            if fnmatch.fnmatch(filename, '*.xml'):
                yield self.process_xml_file(os.path.join(self.directory,filename)), filename

    def write_xml_matches(self,verbosity=0):
        ''' The "real" method:  runs the processing methods and writes the output
            file.  Optionally prints runtime data.
            
            Inputs:  verbosity = 0   Print nothing
                     verbosity = 1   Print filename being processed
                     verbosity = 2   Print matched lines only
                     verbosity = 3   Print both.
        '''
        # Delete the output file if found.
        if os.path.exists(self.output_file):
            os.unlink(self.output_file)

        # Open the output file.
        with codecs.open(self.output_file,'w',encoding='utf-8') as f:
            for match_in_file in self.process_xml_directory():
                start = datetime.datetime.now()
                match_count = 0 
                for match_str in match_in_file[0]:
                    if verbosity in [2,3]:
                        sys.stdout.write(unicode('   %s\n' % match_str.decode('utf-8')))
                        sys.stdout.flush()                          
                    f.write(unicode('%s\n' % match_str.decode('utf-8')))
                    f.flush()  # Keep disk in sync with output.
                    match_count += 1
                if verbosity in [1,3]:
                    sys.stdout.write('(file searched, elapsed time, num_matches, num_unique_pmids): (%s, %s, %s)\n' % 
                                      (match_in_file[1],datetime.datetime.now()-start,match_count))
                    sys.stdout.flush()
            f.close()
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Post-process PMID search into a CSV file.')
    parser.add_argument('-xml', help='Directory of XML files.', required=True)
    parser.add_argument('-tf_list', help='File of transcription factor names and aliases.', required=True)
    parser.add_argument('-match', help='Field to match.', required=True,choices=['xml','abstract', 'title'])
    parser.add_argument('-ofile', help='Output filename.', required=True)
    parser.add_argument('-a', help='Use loose matching for assay names.', action='store_false')

    args = parser.parse_args()

    xml_directory = args.xml
    tf_dict = read_tf_list(args.tf_list)
    match_type = args.match
    output_filename = args.ofile
    strict_assay = args.a
      
    # Adjust assays to remove '+'
    assays = [x.replace('+',' ') for x in assays]
    #for key in tf_dict:
	#	tf_replaced = [x.replace('+',' ') for x in tf_dict[key]]
	#	tf_dict[key]=tf_replaced
		
		
    
    
    procxml = ProcessXmlDirectory(xml_directory,output_filename,match_type,assays,cyts,tf_dict, strict_assay = strict_assay)
    procxml.write_xml_matches(verbosity=1)
    print 'Searching is complete.'
