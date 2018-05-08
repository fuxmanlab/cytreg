# CytReg
This repository contains the python scripts used to retrieve PubMed articles containing at least a combination of a TF, a cytokine, and an assay. Additionally has the cytokine and TF name lists used to query PubMed.

## Description of Files
- nih_eutils_defs.py\
  holds definition of assays and function to create the TF and cytokine dictionaries
- nih_eutils_tf_pmid.py
  downloads PubMed articles XLM files containing at least a TF, cytokine and assay
- post_process_pmids.py
  extracts all TF, cytokine and assays found in a given XML file
  
## Usage

1)  nih_eutils_tf_pmid
    python nih_eutils_tf_pmid.py tf_names.txt cyt_names.txt output_directory
  
2)  post_process_pmids.py
    python post_process_pmids.py -xml ../xml_directory -tf_list tf_names.txt cyt_names.txt -match xml -ofile output_file -a 
