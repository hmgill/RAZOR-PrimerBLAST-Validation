## RAZOR Primer-BLAST Validation Scripts

There is not currently a dedicated API for the Primer-BLAST tool, so we created the scripts below for automatically sending primers to PrimerBLAST and parsing the results.

### Usage 

1. `python3 submit_primerblast_jobs.py`
2. `python3 validate_records_from_primerblast.py primer_jobs_all.json razor_ncbi_primerblast_validation.json`
3. `python3 analyze_validations.py razor_ncbi_primerblast_validation.json`

### File Descriptions
#### razor_db.csv
Copy of all primer pairs for the RAZOR database. Each row represents one primer pair. 

#### submit_primerblast_jobs.py
This script reads primer pairs from `razor_db.csv`, then automatically fills out the Primer-BLAST search form with the appropriate information for the accession, forward primer sequence, and reverse primer sequence. Search parameters are also set according to those specified in the RAZOR manuscript. All jobs submitted to Primer-BLAST are recorded in an output JSON folder. **NOTE: Primer-BLAST retires jobs after 24 hours, so this script will likely need to be run again to create new jobs**

Example Output: 


    "primer_id": "NC_002645.1-135-234-5",
    
    "accession": "NC_002645.1",
    
    "forward_primer": "TGGGTTGCAACAGTTTGGAA",
    
    "reverse_primer": "GGAACCTCGAGTAAGGCGTT",
    
    "product_size": 100,
    
    "left_primer_start": 135,
    
    "right_primer_start": 234,
    
    "job_key": "EhjPIQTRCXkuRwxCASIocHs5OUJWKiJfVw",
    
    "results_url": "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1764392685&job_key=EhjPIQTRCXkuRwxCASIocHs5OUJWKiJfVw",
    
    "submission_time": "2025-11-29 00:04:45",
    
    "status": "submitted"



