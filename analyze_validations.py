import json
from collections import defaultdict
from typing import Dict, List, Tuple

def analyze_primer_validation(json_file_path: str) -> Tuple[int, int, List[Tuple[str, str]]]:
    """
    Analyze primer validation results from a JSON file containing a list of records.
    
    Args:
        json_file_path: Path to the JSON file containing a list of primer records
        
    Returns:
        Tuple containing:
        - Count of all records that passed
        - Count of all records that failed
        - List of tuples (primer_id, validation_status) for all records
    """
    # Track all records
    all_records = []
    passed_count = 0
    failed_count = 0
    
    # Read and parse JSON file
    with open(json_file_path, 'r') as f:
        records = json.load(f)
    
    # Ensure we have a list
    if not isinstance(records, list):
        raise ValueError("JSON file must contain a list of records")
    
    # Process each record
    for record in records:
        primer_id = record.get('primer_id', 'MISSING_ID')
        validation_status = record.get('validation_status')
        
        # Track all occurrences
        all_records.append((primer_id, validation_status))
            
        # Track counts by status (handle both 'fail' and 'failed')
        if validation_status == 'pass':
            passed_count += 1
        elif validation_status in ['fail', 'failed']:
            failed_count += 1
    
    return (passed_count, failed_count, all_records)

def print_results(passed_count: int, failed_count: int, all_records: List[Tuple[str, str]]):
    """Print formatted results of the analysis."""
    print("=" * 60)
    print("PRIMER VALIDATION SUMMARY")
    print("=" * 60)
    print(f"Total validations:  {passed_count + failed_count}")
    print(f"Passed validations: {passed_count}")
    print(f"Failed validations: {failed_count}")
    print("=" * 60)

def save_results_to_file(passed_count: int, failed_count: int,
                         all_records: List[Tuple[str, str]], 
                         output_file: str = 'primer_validation_summary.txt'):
    """Save results to a text file."""
    with open(output_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("PRIMER VALIDATION SUMMARY\n")
        f.write("=" * 60 + "\n")
        f.write(f"Total validations:  {passed_count + failed_count}\n")
        f.write(f"Passed validations: {passed_count}\n")
        f.write(f"Failed validations: {failed_count}\n")
        f.write("=" * 60 + "\n")
    
    print(f"\nResults saved to: {output_file}")

def save_all_records_list(all_records: List[Tuple[str, str]], 
                          output_file: str = 'all_primer_records.txt'):
    """Save all primer IDs with their status to a file."""
    with open(output_file, 'w') as f:
        for primer_id, status in all_records:
            f.write(f"{primer_id}\t{status}\n")
    
    print(f"All records saved to: {output_file}")

if __name__ == "__main__":
    import sys
    
    # Get input file from command line or use default
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = 'primer_data.json'
    
    try:
        # Analyze the data
        passed_count, failed_count, all_records = analyze_primer_validation(input_file)
        
        # Print results to console
        print_results(passed_count, failed_count, all_records)
        
        # Save to file
        save_results_to_file(passed_count, failed_count, all_records)
        
        # Save all records to separate file
        save_all_records_list(all_records)
        
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        print(f"Usage: python {sys.argv[0]} <path_to_json_file>")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON format in file '{input_file}'")
        print(f"Details: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
