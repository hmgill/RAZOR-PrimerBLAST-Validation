## RAZOR Primer-BLAST Validation Scripts

There is not currently a dedicated API for the Primer-BLAST tool, so we created the scripts below for automatically sending primers to PrimerBLAST and parsing the results.

### Usage 

1. `python3 submit_primerblast_jobs.py`
2. `python3 validate_records_from_primerblast.py primer_jobs_all.json razor_ncbi_primerblast_validation.json`
3. `python3 analyze_validations.py razor_ncbi_primerblast_validation.json`

### File Descriptions
#### razor_db.csv
Copy of all primer pairs for the RAZOR database. Each row represents one primer pair. 

#### 1. submit_primerblast_jobs.py
This script reads primer pairs from `razor_db.csv`, then automatically fills out the Primer-BLAST search form with the appropriate information for the accession, forward primer sequence, and reverse primer sequence. Search parameters are also set according to those specified in the RAZOR manuscript. All jobs submitted to Primer-BLAST are recorded in an output JSON folder. 

**note: Primer-BLAST jobs expire after 24 hours, so this script will likely need to be run again to create new jobs**


#### 2. validate_records_from_primerblast.py
The `validate_records_from_primerblast.py` script is used to parse the HTML output from completed jobs. Validation is considered successful if there is a match between the original `razor_db.csv` query pair and a Primer-BLAST result pair. A match occurs when: 
* The forward primer start coordinates for the query / result are identical
* The reverse primer start coordinates for the query / result are identical
* The forward primer sequence for the query / result are identical
* The reverse primer sequence for the query / result are identical
* The PCR product size for the query / result are identical
  
The results are written to a new output JSON file (`razor_validation_results.json`). 

**note: Primer-BLAST jobs take a couple of minutes to run before they are complete. Make sure to give enough time before running this script**

#### 3. analyze_validations.py
This script reads the parsed job results (eg `razor_ncbi_primerblast_validation.json`) and outputs a summary of the numbers of primers that passed / failed validation.
The results are also written to separate text files (`primer_validation_summary.txt`, `all_primer_records.txt`). 

### RESULT SUMMARY 
After running the scripts above, we found that 18826 / 18864 primer pairs passed the Primer-BLAST validation, while 38 did not (`primer_validation_summary.txt`). The errors or failure causes for the rejected primers are shown in the 'failure_cause' field for those entries in `razor_validation_results.json`.




