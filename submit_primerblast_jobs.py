#!/usr/bin/env python3
"""
Script 1 v3: Multithreaded Primer-BLAST Submitter
"""

import pandas as pd
import requests
import json
import time
import re
import logging
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Thread lock for writing to file safely
file_lock = threading.Lock()

class PrimerBlastSubmitter:
    def __init__(self, email: str):
        self.email = email
        self.submit_url = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"
    
    def _determine_product_size_range(self, product_size: int) -> tuple:
        if product_size < 150:
            return (100, 200)
        elif product_size <= 200:
            return (150, 200)
        else:
            return (300, 500)
    
    def submit_job(self, row_data) -> dict:
        """
        Submit a single job. Accepts a dictionary/row to be thread-friendly.
        """
        # Unpack row data
        primer_id = row_data['primer_id']
        accession = row_data['accession']
        forward_primer = row_data['left_primer_seq']
        reverse_primer = row_data['right_primer_seq']
        product_size = int(row_data['product_size'])
        left_primer_start = row_data.get('left_primer_start')
        right_primer_start = row_data.get('right_primer_start')
        organism = "Viruses (taxid:10239)"

        min_size, max_size = self._determine_product_size_range(product_size)
        
        headers = {
            'User-Agent': f'Python-Script/1.0 (mailto:{self.email})',
            'Origin': 'https://www.ncbi.nlm.nih.gov',
            'Referer': 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi'
        }

        params = {
            'CMD': 'request', 
            'INPUT_SEQUENCE': accession,
            'PRIMER_LEFT_INPUT': forward_primer,
            'PRIMER_RIGHT_INPUT': reverse_primer,
            'PRIMER_PRODUCT_MIN': str(min_size),
            'PRIMER_PRODUCT_MAX': str(max_size),
            'PRIMER_MIN_TM': '58',
            'PRIMER_OPT_TM': '60',
            'PRIMER_MAX_TM': '62',
            'PRIMER_MAX_DIFF_TM': '5',
            'ORGANISM': organism,
            'PRIMER_SPECIFICITY_DATABASE': 'nt',
            'TOTAL_PRIMER_SPECIFICITY_MISMATCH': '1',
            'PRIMER_3END_SPECIFICITY_MISMATCH': '1',
            'MISMATCH_REGION_LENGTH': '5',
            'MAX_TARGET_SIZE': str(max_size),
            'SEARCHMODE': '0',
            'SEARCH_SPECIFIC_PRIMER': 'on',
            'ENTREZ_QUERY': accession,
            'NUM_TARGETS': '1000',
            'NUM_TARGETS_WITH_PRIMERS': '1000',
            'SHOW_SVIEWER': 'true',
            'PRIMER_NUM_RETURN': '20',
        }
        
        result_info = {
            'primer_id': primer_id,
            'accession': accession,
            'forward_primer': forward_primer,
            'reverse_primer': reverse_primer,
            'product_size': product_size,
            'left_primer_start': left_primer_start,
            'right_primer_start': right_primer_start,
            'job_key': None,
            'results_url': None,
            'submission_time': time.strftime('%Y-%m-%d %H:%M:%S'),
            'status': 'init'
        }

        try:
            response = requests.post(
                self.submit_url,
                data=params,
                headers=headers,
                allow_redirects=False,
                timeout=60
            )
            response.raise_for_status()
            html_content = response.text
            job_key = None
            results_url = None

            # Strategy 0: META Refresh
            meta_match = re.search(r'content=["\'][^"\']*URL=([^"\']+)["\']', html_content, re.IGNORECASE)
            if meta_match:
                results_url = meta_match.group(1).replace('&amp;', '&')
                key_match = re.search(r'job_key=([A-Za-z0-9_-]+)', results_url)
                if key_match:
                    job_key = key_match.group(1)

            # Strategy 1: Hidden input
            if not job_key:
                input_match = re.search(r'name="job_key"\s+value="([^"]+)"', html_content)
                if not input_match:
                     input_match = re.search(r'value="([^"]+)"\s+name="job_key"', html_content)
                if input_match:
                    job_key = input_match.group(1)

            # Strategy 2: Text display
            if not job_key:
                job_id_match = re.search(r'JOB\s+ID\s*:\s*([A-Za-z0-9_-]+)', html_content, re.IGNORECASE)
                if job_id_match:
                    job_key = job_id_match.group(1)
            
            if job_key and not results_url:
                results_url = f"{self.submit_url}?job_key={job_key}"

            # Update result info
            result_info['job_key'] = job_key
            result_info['results_url'] = results_url
            result_info['status'] = 'submitted' if job_key else 'failed_extraction'

            # REMOVED: The block that saved the HTML file on failure
            
        except Exception as e:
            logger.error(f"Error submitting {primer_id}: {e}")
            result_info['status'] = 'failed'
            result_info['error'] = str(e)
            
        return result_info

def submit_primer_jobs(csv_path: str, email: str, output_path: str = None, start_line: int = 0, end_line: int = None):
    submitter = PrimerBlastSubmitter(email=email)
    df = pd.read_csv(csv_path)
    
    total_rows = len(df)
    actual_end = end_line if end_line is not None else total_rows
    df_subset = df.iloc[start_line:actual_end].reset_index(drop=True)
    
    if output_path is None:
        range_str = f"all" if end_line is None else f"{start_line}_to_{end_line}"
        output_path = f"primer_jobs_{range_str}.json"
    
    print(f"Processing rows {start_line} to {actual_end-1} ({len(df_subset)} jobs)")
    print(f"Using Multithreading (Max 4 workers to respect NCBI limits)...")
    
    completed_jobs = []
    
    # Initialize file with empty list if it doesn't exist
    with open(output_path, 'w') as f:
        json.dump([], f)

    # Use ThreadPoolExecutor
    # max_workers=4 is safe. NCBI allows ~3 requests/sec. 
    # With network overhead, 4 workers usually averages out to safe limits.
    with ThreadPoolExecutor(max_workers=8) as executor:
        # Create a list of futures
        futures = []
        for _, row in df_subset.iterrows():
            futures.append(executor.submit(submitter.submit_job, row))
        
        # Process as they complete
        count = 0
        for future in as_completed(futures):
            result = future.result()
            count += 1
            
            # Print status
            symbol = "✓" if result['job_key'] else "✗"
            print(f"[{count}/{len(df_subset)}] {symbol} {result['primer_id']} -> {result['job_key'] or 'FAILED'}")
            
            # Thread-safe append to memory list
            completed_jobs.append(result)
            
            # Periodically save to disk (every 10 jobs) to prevent data loss on crash
            if count % 10 == 0 or count == len(df_subset):
                with file_lock:
                    with open(output_path, 'w') as f:
                        json.dump(completed_jobs, f, indent=2)

    print(f"\nAll jobs processed. Results saved to: {output_path}")

if __name__ == "__main__":
    # Configuration
    EMAIL = "" ## Your E-Mail Here
    INPUT_CSV = "razor_db.csv"
    
    # CLI Arguments logic preserved
    START_LINE = 0
    END_LINE = None
    
    if len(sys.argv) > 1:
        START_LINE = int(sys.argv[1])
    if len(sys.argv) > 2:
        END_LINE = int(sys.argv[2])
    
    submit_primer_jobs(INPUT_CSV, EMAIL, None, START_LINE, END_LINE)
