#!/usr/bin/env python3
"""
Script 2: Validate completed Primer-BLAST jobs (Multithreaded Version)
"""

import json
import requests
import re
import time
import logging
from typing import Dict, Optional, List
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
import sys

# Configure logging
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class PrimerBlastValidator:
    """Validate completed Primer-BLAST jobs with thread-safe rate limiting"""
    
    def __init__(self, max_workers: int = 8):
        """
        Initialize validator with thread-safe rate limiting
        
        Args:
            max_workers: Maximum number of concurrent threads (default: 3)
                        NCBI allows ~3 requests/second, so 3 workers is safe
        """
        self.last_request_time = 0
        self.rate_limit_lock = Lock()
        self.max_workers = max_workers
    
    def _rate_limit(self):
        """Thread-safe rate limiting to respect NCBI's 3 requests/second limit"""
        with self.rate_limit_lock:
            current_time = time.time()
            time_since_last = current_time - self.last_request_time
            if time_since_last < 0.34:  # ~3 requests per second
                time.sleep(0.34 - time_since_last)
            self.last_request_time = time.time()
    
    def fetch_results(self, results_url: str) -> Optional[str]:
        """Fetch HTML results from Primer-BLAST"""
        self._rate_limit()
        
        try:
            response = requests.get(results_url, timeout=60)
            response.raise_for_status()
            return response.text
        except Exception as e:
            logger.error(f"Error fetching results: {e}")
            return None
    
    def extract_all_primer_pairs(self, html_content: str) -> List[Dict]:
        """
        Extract all primer pairs from HTML when multiple pairs are present
        
        Returns:
            List of dicts, each containing primer info for one pair
        """
        all_pairs = []
        
        # Split HTML into primer pair sections
        # Look for <h3>Primer pair N</h3> headers
        pair_sections = re.split(r'<h3[^>]*>Primer pair \d+</h3>', html_content, flags=re.IGNORECASE)
        
        if len(pair_sections) <= 1:
            # No explicit primer pair headers, treat entire HTML as one pair
            pair_sections = [html_content]
        else:
            # Remove the first element (everything before first primer pair)
            pair_sections = pair_sections[1:]
        
        logger.debug(f"Found {len(pair_sections)} primer pair sections")
        
        for idx, section in enumerate(pair_sections, 1):
            pair_data = self._extract_single_pair(section, pair_index=idx)
            if pair_data['forward_start'] is not None:  # Only add if we found data
                all_pairs.append(pair_data)
        
        return all_pairs
    
    def _extract_single_pair(self, html_content: str, pair_index: int = None) -> Dict:
        """
        Extract primer coordinates and product length from a single primer pair section
        
        Looks for table rows like:
        <tr><th>Forward primer</th><td>TTTTTGCTATGGCCGGCATC</td><td>Plus</td><td>20</td><td>78</td><td>97</td><td>59.54</td><td>50.00</td><td>6.00</td><td>4.00</td></tr>
        <tr><th>Reverse primer</th><td>TGATTGCCATGCAGACCCAT</td><td>Minus</td><td>20</td><td>498</td><td>479</td><td>59.54</td><td>50.00</td><td>6.00</td><td>4.00</td></tr>
        <tr><th>Product length</th><td colspan="9">421</td></tr>
        
        Columns: Sequence, Strand, Length, Start, Stop, Tm, GC%, Self-complementarity, Self-3' complementarity
        
        Returns:
            Dict with forward_start, forward_end, reverse_start, reverse_end, product_length, sequences
        """
        result = {
            'forward_start': None,
            'forward_end': None,
            'forward_sequence': None,
            'reverse_start': None,
            'reverse_end': None,
            'reverse_sequence': None,
            'product_length': None,
            'pair_index': pair_index
        }
        
        # Parse Forward primer row
        forward_pattern = r'<tr><th>Forward\s+primer</th><td>([ATCGN]+)</td><td>([^<]+)</td><td>(\d+)</td><td>(\d+)</td><td>(\d+)</td>'
        forward_match = re.search(forward_pattern, html_content, re.IGNORECASE)
        
        if forward_match:
            result['forward_sequence'] = forward_match.group(1)
            strand = forward_match.group(2)
            length = int(forward_match.group(3))
            start = int(forward_match.group(4))
            stop = int(forward_match.group(5))
            
            result['forward_start'] = start
            result['forward_end'] = stop
            
            logger.debug(f"Extracted forward (pair {pair_index}): seq={result['forward_sequence']}, start={start}, end={stop}, strand={strand}")
        else:
            logger.warning(f"Could not find Forward primer row in HTML (pair {pair_index})")
        
        # Parse Reverse primer row
        reverse_pattern = r'<tr><th>Reverse\s+primer</th><td>([ATCGN]+)</td><td>([^<]+)</td><td>(\d+)</td><td>(\d+)</td><td>(\d+)</td>'
        reverse_match = re.search(reverse_pattern, html_content, re.IGNORECASE)
        
        if reverse_match:
            result['reverse_sequence'] = reverse_match.group(1)
            strand = reverse_match.group(2)
            length = int(reverse_match.group(3))
            start = int(reverse_match.group(4))
            stop = int(reverse_match.group(5))
            
            result['reverse_start'] = start
            result['reverse_end'] = stop
            
            logger.debug(f"Extracted reverse (pair {pair_index}): seq={result['reverse_sequence']}, start={start}, end={stop}, strand={strand}")
        else:
            logger.warning(f"Could not find Reverse primer row in HTML (pair {pair_index})")
        
        # Parse Product length
        product_pattern = r'<tr><th>Product\s+length</th><td[^>]*>(\d+)</td></tr>'
        product_match = re.search(product_pattern, html_content, re.IGNORECASE)
        
        if product_match:
            result['product_length'] = int(product_match.group(1))
            logger.debug(f"Extracted product length (pair {pair_index}): {result['product_length']}")
        else:
            logger.warning(f"Could not find Product length row in HTML (pair {pair_index})")
        
        return result
    
    def extract_coordinates(self, html_content: str) -> Dict:
        """
        Extract primer coordinates, scanning through all primer pairs if multiple exist.
        For backward compatibility, returns the first pair found.
        
        Returns:
            Dict with forward_start, forward_end, reverse_start, reverse_end, product_length, sequences
        """
        all_pairs = self.extract_all_primer_pairs(html_content)
        
        if not all_pairs:
            return {
                'forward_start': None,
                'forward_end': None,
                'forward_sequence': None,
                'reverse_start': None,
                'reverse_end': None,
                'reverse_sequence': None,
                'product_length': None,
                'pair_index': None
            }
        
        # Return first pair for backward compatibility
        return all_pairs[0]
    
    def find_matching_pair(self, all_pairs: List[Dict], job_info: Dict) -> tuple:
        """
        Find the primer pair that matches the expected coordinates
        
        Args:
            all_pairs: List of all extracted primer pairs
            job_info: Job information with expected values
        
        Returns:
            Tuple of (matched_pair_dict, match_found_boolean)
        """
        expected_fwd = job_info.get('left_primer_start')
        expected_rev = job_info.get('right_primer_start')
        expected_size = job_info.get('product_size')
        
        for pair in all_pairs:
            # Check if this pair matches
            fwd_match = (expected_fwd is None or pair['forward_start'] is None or 
                        int(expected_fwd) == int(pair['forward_start']))
            rev_match = (expected_rev is None or pair['reverse_start'] is None or 
                        int(expected_rev) == int(pair['reverse_start']))
            size_match = (expected_size is None or pair['product_length'] is None or 
                         int(expected_size) == int(pair['product_length']))
            
            if fwd_match and rev_match and size_match:
                logger.debug(f"Found matching pair at index {pair['pair_index']}")
                return pair, True
        
        # No match found, return first pair
        return all_pairs[0] if all_pairs else None, False
    
    def validate_job(self, job_info: Dict, job_index: int = None, total_jobs: int = None) -> Dict:
        """
        Validate a single job by fetching results and comparing coordinates
        
        Args:
            job_info: Job information from submit_jobs.py
            job_index: Index of this job in the batch (for progress display)
            total_jobs: Total number of jobs being processed
        
        Returns:
            Updated job_info with validation results
        """
        if not job_info.get('results_url'):
            job_info['validation_status'] = 'no_url'
            if job_index is not None:
                print(f"[{job_index}/{total_jobs}] {job_info['primer_id']} - ✗ No URL")
            return job_info
        
        if job_index is not None:
            print(f"[{job_index}/{total_jobs}] {job_info['primer_id']} - Fetching...", flush=True)
        
        html = self.fetch_results(job_info['results_url'])
        
        if not html:
            if job_index is not None:
                print(f"[{job_index}/{total_jobs}] {job_info['primer_id']} - ✗ Fetch failed")
            job_info['validation_status'] = 'fetch_failed'
            return job_info
        
        # Check if job is complete
        if 'processing' in html.lower() or 'queued' in html.lower():
            if job_index is not None:
                print(f"[{job_index}/{total_jobs}] {job_info['primer_id']} - ⏳ Still processing")
            job_info['validation_status'] = 'processing'
            return job_info
        
        if 'No primer pairs found' in html:
            if job_index is not None:
                print(f"[{job_index}/{total_jobs}] {job_info['primer_id']} - ✗ No primers found")
            job_info['validation_status'] = 'no_primers'
            return job_info
        
        # Extract all primer pairs
        all_pairs = self.extract_all_primer_pairs(html)
        
        if not all_pairs:
            if job_index is not None:
                print(f"[{job_index}/{total_jobs}] {job_info['primer_id']} - ✗ Could not extract primers")
            job_info['validation_status'] = 'extraction_failed'
            return job_info
        
        # Record total number of pairs found
        job_info['total_primer_pairs_found'] = len(all_pairs)
        
        # Find matching pair or use first one
        matched_pair, match_found = self.find_matching_pair(all_pairs, job_info)
        
        if matched_pair is None:
            if job_index is not None:
                print(f"[{job_index}/{total_jobs}] {job_info['primer_id']} - ✗ No valid primer data")
            job_info['validation_status'] = 'no_primer_data'
            return job_info
        
        # Store which pair was used for validation
        job_info['matched_pair_index'] = matched_pair.get('pair_index')
        job_info['used_matching_algorithm'] = match_found
        
        # Store extracted data
        job_info['extracted_forward_sequence'] = matched_pair['forward_sequence']
        job_info['extracted_forward_start'] = matched_pair['forward_start']
        job_info['extracted_forward_end'] = matched_pair['forward_end']
        job_info['extracted_reverse_sequence'] = matched_pair['reverse_sequence']
        job_info['extracted_reverse_start'] = matched_pair['reverse_start']
        job_info['extracted_reverse_end'] = matched_pair['reverse_end']
        job_info['extracted_product_length'] = matched_pair['product_length']
        
        # Check if sequences match
        sequence_match = True
        if matched_pair['forward_sequence'] and job_info.get('forward_primer'):
            if matched_pair['forward_sequence'].upper() != job_info['forward_primer'].upper():
                sequence_match = False
                logger.warning(f"Forward sequence mismatch: expected {job_info['forward_primer']}, got {matched_pair['forward_sequence']}")
        
        if matched_pair['reverse_sequence'] and job_info.get('reverse_primer'):
            if matched_pair['reverse_sequence'].upper() != job_info['reverse_primer'].upper():
                sequence_match = False
                logger.warning(f"Reverse sequence mismatch: expected {job_info['reverse_primer']}, got {matched_pair['reverse_sequence']}")
        
        job_info['sequences_match'] = sequence_match
        
        # Validate coordinates
        validation = {
            'forward_start_match': False,
            'reverse_start_match': False,
            'product_size_match': False,
            'sequences_match': sequence_match
        }
        
        expected_fwd = job_info.get('left_primer_start')
        expected_rev = job_info.get('right_primer_start')
        expected_size = job_info.get('product_size')
        
        if expected_fwd is not None and matched_pair['forward_start'] is not None:
            validation['forward_start_match'] = (
                int(expected_fwd) == int(matched_pair['forward_start'])
            )
        
        if expected_rev is not None and matched_pair['reverse_start'] is not None:
            validation['reverse_start_match'] = (
                int(expected_rev) == int(matched_pair['reverse_start'])
            )
        
        if expected_size is not None and matched_pair['product_length'] is not None:
            validation['product_size_match'] = (
                int(expected_size) == int(matched_pair['product_length'])
            )
        
        job_info['validation'] = validation
        job_info['validation_passed'] = (
            validation['forward_start_match'] and
            validation['reverse_start_match'] and
            validation['product_size_match'] and
            validation['sequences_match']
        )
        
        # Prepare status message
        status_msg = ""
        if len(all_pairs) > 1:
            status_msg = f" ({len(all_pairs)} pairs, using #{matched_pair.get('pair_index', 1)})"
        
        if job_info['validation_passed']:
            if job_index is not None:
                print(f"[{job_index}/{total_jobs}] {job_info['primer_id']}{status_msg} - ✓ PASS")
            job_info['validation_status'] = 'pass'
        else:
            failed = []
            if not validation['sequences_match']:
                failed.append("SEQ_MISMATCH")
            if not validation['forward_start_match']:
                failed.append(f"F_start:{expected_fwd}≠{matched_pair['forward_start']}")
            if not validation['reverse_start_match']:
                failed.append(f"R_start:{expected_rev}≠{matched_pair['reverse_start']}")
            if not validation['product_size_match']:
                failed.append(f"Size:{expected_size}≠{matched_pair['product_length']}")
            if job_index is not None:
                print(f"[{job_index}/{total_jobs}] {job_info['primer_id']}{status_msg} - ✗ FAIL ({', '.join(failed)})")
            job_info['validation_status'] = 'fail'
        
        return job_info


def validate_all_jobs(
    jobs_path: str = "primer_jobs_all.json",
    output_path: str = None,
    verbose: bool = False,
    start_line: int = 0,
    end_line: int = None,
    max_workers: int = 8,
    checkpoint_interval: int = 1000
):
    """
    Validate all submitted jobs using multithreading
    
    Args:
        jobs_path: Path to JSON file from submit_jobs.py
        output_path: Where to save validation results (if None, auto-generated based on line range)
        verbose: Show detailed extraction information
        start_line: Start validating from this job number (0-indexed, default: 0)
        end_line: End validating at this job number (exclusive, default: None = validate to end)
        max_workers: Maximum number of concurrent threads (default: 3)
        checkpoint_interval: Save results every N completed jobs (default: 10)
    """
    # Load jobs
    with open(jobs_path) as f:
        jobs = json.load(f)
    
    # Determine the range to process
    total_jobs = len(jobs)
    actual_end = end_line if end_line is not None else total_jobs
    
    # Validate range
    if start_line < 0:
        print(f"Error: start_line ({start_line}) cannot be negative")
        return []
    if start_line >= total_jobs:
        print(f"Error: start_line ({start_line}) is beyond the jobs list length ({total_jobs} jobs)")
        return []
    if actual_end > total_jobs:
        print(f"Warning: end_line ({end_line}) is beyond jobs list length ({total_jobs}), using {total_jobs}")
        actual_end = total_jobs
    if start_line >= actual_end:
        print(f"Error: start_line ({start_line}) must be less than end_line ({actual_end})")
        return []
    
    # Slice the jobs list
    jobs_subset = jobs[start_line:actual_end]
    
    # Auto-generate output filename if not provided
    if output_path is None:
        if end_line is None:
            output_path = f"validation_results_{start_line}_to_end.json"
        else:
            output_path = f"validation_results_{start_line}_to_{end_line}.json"
    
    print(f"Validating jobs {start_line} to {actual_end-1} ({len(jobs_subset)} total jobs)")
    print(f"Using {max_workers} worker threads")
    print(f"Checkpointing every {checkpoint_interval} completed jobs")
    if verbose:
        print("(Verbose mode: showing extraction details)")
    print("="*80)
    
    validator = PrimerBlastValidator(max_workers=max_workers)
    validated_jobs = [None] * len(jobs_subset)  # Pre-allocate to maintain order
    completed_count = 0
    save_lock = Lock()
    
    def process_job(idx: int, job: Dict) -> tuple:
        """Process a single job and return its index and result"""
        actual_job_num = idx + start_line
        
        if not job.get('job_key'):
            print(f"[{idx + 1}/{len(jobs_subset)} (job {actual_job_num})] {job['primer_id']} - ✗ No job key")
            job['validation_status'] = 'no_job_key'
            return idx, job
        
        validated = validator.validate_job(
            job,
            job_index=idx + 1,
            total_jobs=len(jobs_subset)
        )
        
        # Show extraction details in verbose mode
        if verbose and validated.get('extracted_forward_start') is not None:
            print(f"  Expected: F={job.get('left_primer_start')}, R={job.get('right_primer_start')}, Size={job.get('product_size')}")
            print(f"  Extracted: F={validated.get('extracted_forward_start')}, R={validated.get('extracted_reverse_start')}, Size={validated.get('extracted_product_length')}")
            if validated.get('extracted_forward_sequence'):
                print(f"  Forward seq: {validated['extracted_forward_sequence']}")
            if validated.get('extracted_reverse_sequence'):
                print(f"  Reverse seq: {validated['extracted_reverse_sequence']}")
        
        return idx, validated
    
    # Use ThreadPoolExecutor for concurrent validation
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_idx = {
            executor.submit(process_job, idx, job): idx
            for idx, job in enumerate(jobs_subset)
        }
        
        # Process completed jobs as they finish
        for future in as_completed(future_to_idx):
            idx, validated_job = future.result()
            validated_jobs[idx] = validated_job
            completed_count += 1
            
            # Checkpoint: save progress periodically
            if completed_count % checkpoint_interval == 0:
                with save_lock:
                    # Filter out None values (not yet completed)
                    jobs_to_save = [j for j in validated_jobs if j is not None]
                    with open(output_path, 'w') as f:
                        json.dump(jobs_to_save, f, indent=2)
                    print(f"\n[CHECKPOINT] Saved {completed_count}/{len(jobs_subset)} completed jobs to {output_path}\n")
    
    # Final save with all completed jobs
    with open(output_path, 'w') as f:
        json.dump(validated_jobs, f, indent=2)
    
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    
    # Show range info
    if start_line > 0 or end_line is not None:
        print(f"\nProcessed jobs {start_line} to {actual_end-1}")
    
    # Count statuses
    statuses = {}
    for job in validated_jobs:
        status = job.get('validation_status', 'unknown')
        statuses[status] = statuses.get(status, 0) + 1
    
    print(f"\nTotal jobs: {len(validated_jobs)}")
    print(f"\nStatus breakdown:")
    for status, count in sorted(statuses.items()):
        print(f"  {status}: {count}")
    
    if 'pass' in statuses:
        print(f"\n✓ Validation passed: {statuses['pass']}/{len(validated_jobs)}")
    
    if 'fail' in statuses:
        print(f"✗ Validation failed: {statuses['fail']}/{len(validated_jobs)}")
    
    if 'processing' in statuses:
        print(f"⏳ Still processing: {statuses['processing']}/{len(validated_jobs)}")
        print("\n  Wait longer and run script again")
    
    # Statistics about multiple primer pairs
    jobs_with_multiple_pairs = sum(1 for job in validated_jobs 
                                   if job.get('total_primer_pairs_found', 0) > 1)
    jobs_using_matching = sum(1 for job in validated_jobs 
                             if job.get('used_matching_algorithm', False))
    
    if jobs_with_multiple_pairs > 0:
        print(f"\n📊 Multiple primer pairs detected in {jobs_with_multiple_pairs} jobs")
        print(f"   Matching algorithm used successfully: {jobs_using_matching} jobs")
        
        # Show distribution of pair counts
        pair_counts = {}
        for job in validated_jobs:
            count = job.get('total_primer_pairs_found', 0)
            if count > 0:
                pair_counts[count] = pair_counts.get(count, 0) + 1
        
        if len(pair_counts) > 1:
            print(f"\n   Distribution of primer pairs found:")
            for pair_count in sorted(pair_counts.keys()):
                print(f"     {pair_count} pair(s): {pair_counts[pair_count]} jobs")
    
    print(f"\nResults saved to: {output_path}")
    
    return validated_jobs


if __name__ == "__main__":
    # Default configuration
    jobs_path = "primer_jobs_all.json"
    verbose = '--verbose' in sys.argv or '-v' in sys.argv
    start_line = 0
    end_line = None
    max_workers = 8
    checkpoint_interval = 1000
    
    # Parse command-line arguments
    # Usage: python validate_results.py [start_line] [end_line] [-v] [--workers N] [--checkpoint N]
    # Example: python validate_results.py 100 200 -v --workers 5 --checkpoint 20
    args = [arg for arg in sys.argv[1:] if arg not in ['--verbose', '-v']]
    
    # Parse --workers flag
    if '--workers' in sys.argv:
        workers_idx = sys.argv.index('--workers')
        if workers_idx + 1 < len(sys.argv):
            try:
                max_workers = int(sys.argv[workers_idx + 1])
                print(f"Using {max_workers} worker threads")
                # Remove from args list
                args = [arg for i, arg in enumerate(args) if arg != '--workers' and (i == 0 or args[i-1] != '--workers')]
            except ValueError:
                print(f"Error: Invalid --workers value. Using default (3)")
    
    # Parse --checkpoint flag
    if '--checkpoint' in sys.argv:
        checkpoint_idx = sys.argv.index('--checkpoint')
        if checkpoint_idx + 1 < len(sys.argv):
            try:
                checkpoint_interval = int(sys.argv[checkpoint_idx + 1])
                print(f"Checkpointing every {checkpoint_interval} jobs")
                # Remove from args list
                args = [arg for i, arg in enumerate(args) if arg != '--checkpoint' and (i == 0 or args[i-1] != '--checkpoint')]
            except ValueError:
                print(f"Error: Invalid --checkpoint value. Using default (10)")
    
    if len(args) > 0:
        try:
            start_line = int(args[0])
            print(f"Command-line: starting from job {start_line}")
        except ValueError:
            print(f"Error: Invalid start_line value '{args[0]}'. Using default (0)")
    
    if len(args) > 1:
        try:
            end_line = int(args[1])
            print(f"Command-line: ending at job {end_line-1} (exclusive end: {end_line})")
        except ValueError:
            print(f"Error: Invalid end_line value '{args[1]}'. Validating to end")
            end_line = None
    
    results = validate_all_jobs(
        jobs_path=jobs_path,
        output_path=None,
        verbose=verbose,
        start_line=start_line,
        end_line=end_line,
        max_workers=max_workers,
        checkpoint_interval=checkpoint_interval
    )
