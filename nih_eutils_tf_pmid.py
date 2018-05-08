# -*- coding: utf-8 -*-
 
 
from bs4 import BeautifulSoup

import requests

import time
import datetime
import csv
import os
import sys

reload(sys)
sys.setdefaultencoding('utf-8')





# Another stab at a better organized URL retrieval scheme for NIH Eutils

# Load cyts and assays
from nih_eutils_defs_alvaro import *

######### These values are registered with the NIH            ######### 
EUTILS_TOOL = 'nih_eutils_tf_pmid.py'
EUTILS_EMAIL = 'scp@bu.edu'
######### see: https://www.ncbi.nlm.nih.gov/books/NBK25497/   ######### 


def chunks(l, n):
    '''Yield successive n-sized chunks from l.
       from: http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks '''
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


class PmidQuery(object):
    ''' This is a class that takes in a query set of transcription factors,
        cytokines, and assay names.  It provides a method to retrieve
        any matching PMIDs from the query set.  Errors will throw an exception. '''
	

    # URL used to retrieve document abstracts
    NIH_FETCH_URL='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    # URL to perform searches
    NIH_SEARCH_URL='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    # Maximum number of records to return for a query.  This should be
    # much bigger than the anticipated results. Make sure this is tested
    # to double check that Eutils won't barf on it.   
    RETMAX = 8000
    # Size of search terms to be combined into a single query chunk.  Small enough
    # to never cause Eutils problems, large enough to reduce the number of queries.
    CHUNK_LENGTH = 800
    # Set timeouts for the URLs
    FETCH_TIMEOUT = 45
    # These are large queries so they need a big timeout.
    SEARCH_TIMEOUT = 10


    def __init__(self,tf_dict, cyt_dict, assay_list, output_dir):
        ''' Simple constructor '''
        self.tf_dict = tf_dict
        self.cyt_dict = cyt_dict
        self.assay_list = assay_list
        # Create a private session for url fetching
        self.session = requests.Session()
        # To be filled in...
        self.webenv = None
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
    
        
    def fetch_pmid_records(self):
        ''' Externally called method to do a search and then
            retrieve any matching PMIDs '''
        # 1st: run a query over the cytokines and get a final query key.
        cyt_query_key = ''#self._search_cytokines()
        # 2nd: run a query over assays.
        assay_query_key = self._search_assays()
        # 3rd: run a query over transcription factors
        tf_query_key = ''#self._search_transcription_factors()
        print '%s    %s  %s   %s' % (cyt_query_key,assay_query_key,tf_query_key,self.webenv)
        # Now combine the query keys with an AND and get the final record count.
        key_query = ['#%s' % x for x in (cyt_query_key, assay_query_key, tf_query_key)]
        combo = self._search_list_by_chunks(key_query,boolean_op='AND')
        combo_query_key = combo['query_key']
        combo_count = combo['count']
        print 'Combo key: %s    count: %s' % (combo_query_key,combo_count)
        self._fetch_pmids(combo_query_key,combo_count)


    def _search_cytokines(self):
        print 'Building cytokine query set'
        return self._search_dict_by_chunks(self.cyt_dict)['query_key']
        
    def _search_transcription_factors(self):
        print 'Building TF query set'
        return self._search_dict_by_chunks(self.tf_dict)['query_key']
        
    def _search_assays(self):
        print 'Building assay query set'
        return self._search_dict_by_chunks({'assays':self.assay_list})['query_key']     

    def _search_dict_by_chunks(self,the_dict,boolean_op='OR'):
        search_list=[]
        for c in the_dict:
          search_list += the_dict[c]
        return self._search_list_by_chunks(search_list,boolean_op)
    
    def _search_list_by_chunks(self,search_list,boolean_op='OR'):
        # Loop over this in chunks of self.CHUNK_LENGTH, query, and add them to prev queries
        # as we go via OR.
        bad_chunks=[]
        last_query_key=None
        count = 0
        num_chunks = len(search_list) / self.CHUNK_LENGTH
        if len(search_list) % self.CHUNK_LENGTH != 0:
            num_chunks += 1
        
        for i,chunk in enumerate(chunks(search_list,self.CHUNK_LENGTH)):
            try:
                op = ' %s ' % boolean_op
                search_query = op.join(chunk)
                params = {'db':'pubmed', 'usehistory':'y','retmax':0}
                if self.webenv:
                    params['WebEnv'] = self.webenv
                params['term'] = search_query
                if last_query_key:
                    # Use the last result and OR it to this one.
                    params['term'] = '#%s OR %s' % (last_query_key,params['term'])
                url_text = self._retrieve_url(self.NIH_SEARCH_URL,params,timeout=self.SEARCH_TIMEOUT)
                soup = BeautifulSoup(url_text, "lxml")
                last_query_key = soup.find('querykey').get_text()
                if not self.webenv:
                    self.webenv = soup.find('webenv').get_text()
                count = int(soup.find('count').get_text())
                print 'Query chunk %s of %s.    Count of results: %s' % (i+1,num_chunks,count)
                sys.stdout.flush()
            except Exception,e:
                print 'Exception! %s' % e
                bad_chunks.append( (i,e) ) 
        if len(bad_chunks) > 0:
            print 'Query had %s chunks fail to communicate with the server.' % (len(bad_chunks))
        return {'query_key':last_query_key,'count':count}
        
       
    def _fetch_pmids(self,query_key,query_count):
        ''' Internal method.  Using a query_key and the stored webenv,
            it will download full records into XML format.
        '''
        # This will fetch records in pages of self.RETMAX.       
        num_pages = query_count / self.RETMAX
        if query_count % self.RETMAX != 0:
            num_pages += 1
        bad_pages=[]
        
        for i in xrange(num_pages):
            try:
                # The retstart parameter is used to get self.RETMAX records at a time.
                #params = {'db':'pubmed', 'usehistory':'y','retmax':self.RETMAX,
                #         'retmode':'xml','rettype':'abstract','retstart':i * self.RETMAX,
                #         'WebEnv':self.webenv,'query_key':query_key}
                params = {'db':'pubmed', 'usehistory':'y','retmax':self.RETMAX,
                         'retmode':'xml','retstart':i * self.RETMAX,
                         'WebEnv':self.webenv,'query_key':query_key}
                url_text = self._retrieve_url(self.NIH_FETCH_URL,params,timeout=self.FETCH_TIMEOUT)
    
                #####################################
                ### Write files to XML            ###
                fname = os.path.join(self.output_dir, 'page_%s.xml' % i)
                with open(fname,'w') as f:
                    f.write(url_text)
                    f.close()
                #####################################
                print 'Page %s of %s retrieved.' % (i+1,num_pages)
                sys.stdout.flush()
            except Exception,e:
                print 'Exception! %s' % e
                bad_pages.append( (i,e) ) 
        if len(bad_pages) > 0:
            print 'Query had %s pages fail to communicate with the server.' % (len(bad_pages))
        
        return 
    
    
    def _get_assays_query(assays):
      ''' Transform the list of assays into a Boolean search query '''
      return ' OR '.join(['(%s)' % (assay.replace('+',' AND ')) for assay in assays])
    
    def _retrieve_url(self,url,params,timeout=1.0,max_tries=10):
        ''' Calls a url with a timeout.  Will retry if there
            is an error.  '''
        counter = 0
        error = None 
        
        # Add the EUtils email and tool parameters
        params['email'] = EUTILS_EMAIL
        params['tool'] = EUTILS_TOOL
                
        while counter < max_tries:
            if counter > 0:
                msg = 'Retry #%s for URL: %s.' % (counter,url)
                if error:
                    msg += '\n EXCEPTION: %s' % error
                print msg
                
            # Pause to avoid hammering the NIH server!  3x/sec is the max allowed rate.
            time.sleep(0.4)
            try:
                url_response = self.session.post(url,params=params, timeout=timeout) 
                # If everything is cool then return the text.
                if url_response.status_code == requests.codes.ok:
                    url_response.encoding = 'utf-8'    #### NEW LINE TO ADD
                    return url_response.text

            except Exception, e:
                error = e
            counter += 1
            timeout = timeout * 1.333  # Increase the timeout
        # Got to here?  There was an error.  Re-throw the exception.
        raise e


if __name__=='__main__':
    if len(sys.argv) != 4:
        print 'Usage:  python nih_eutils_tf_pmid.py tf_filename cyt_filename output_dir'
        exit(1)
    tf_dict = read_tf_list(sys.argv[1])
    cyts = read_tf_list(sys.argv[2])
    output_dir = sys.argv[3]
    pmq = PmidQuery(tf_dict, cyts, assays,output_dir) 
    start = datetime.datetime.now()
    pmq.fetch_pmid_records()
    end = datetime.datetime.now()
    print 'ELAPSED: %s' % (end-start)
