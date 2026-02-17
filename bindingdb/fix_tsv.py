#!/usr/bin/env python3

import sys
import argparse
import pathlib
import os
import duckdb
from tqdm.auto import tqdm

def fix_tsv(infile_path, tempfile_path, total_size):
    """
    Reads a potentially malformed TSV and writes a normalized version
    where every line has the same number of columns as the header.
    """
    num_columns = 0
    # Using errors='replace' to handle potential encoding issues common in large datasets
    with open(infile_path, 'r', encoding='utf-8', errors='replace') as ifs, \
         open(tempfile_path, 'w', encoding='utf-8') as ofs:
        
        pbar = tqdm(total=total_size, unit='B', unit_scale=True, desc="Fixing TSV")
        
        # Process header
        header = ifs.readline()
        if not header:
            pbar.close()
            return 0
        
        pbar.update(len(header)) # Approximation of bytes
        header_toks = header.rstrip('\r\n').split('\t')
        num_columns = len(header_toks)
        ofs.write('\t'.join(header_toks) + '\n')
        
        # Process data rows
        for line in ifs:
            pbar.update(len(line))
            
            toks = line.rstrip('\r\n').split('\t')
            num_toks = len(toks)
            
            if num_toks < num_columns:
                # Pad missing columns
                toks.extend([''] * (num_columns - num_toks))
            elif num_toks > num_columns:
                # Truncate extra columns
                toks = toks[:num_columns]
            
            ofs.write('\t'.join(toks) + '\n')
            
        pbar.close()
    return num_columns

def main():
    parser = argparse.ArgumentParser(description="Normalize BindingDB TSV and load into DuckDB.")
    parser.add_argument("infile", type=pathlib.Path, help="Input TSV file")
    parser.add_argument("outfile", type=pathlib.Path, help="Output DuckDB file (.ddb)")
    parser.add_argument("--table", default="bindingdb", help="Table name in DuckDB (default: bindingdb)")
    parser.add_argument("--tempfile", type=pathlib.Path, help="Path for temporary fixed TSV (optional)")
    parser.add_argument("--keep-temp", action="store_true", help="Do not delete the temporary file after import")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()

    if args.infile == args.outfile:
        print("Error: Input and output filenames must differ.")
        sys.exit(1)

    if args.outfile.exists():
        print(f"Error: Output file '{args.outfile}' already exists. Please remove or rename it.")
        sys.exit(1)

    # Use provided tempfile or create one based on output name
    temp_path = args.tempfile or args.outfile.with_suffix('.fixed.tsv')
    
    if temp_path.exists():
        print(f"Warning: Temporary file '{temp_path}' already exists and will be overwritten.")

    total_size = args.infile.stat().st_size
    
    try:
        print(f"Fixing '{args.infile}'...")
        num_cols = fix_tsv(args.infile, temp_path, total_size)
        
        if num_cols == 0:
            print("Error: Input file appears to be empty.")
            sys.exit(1)
            
        print(f"Loading into DuckDB (table: '{args.table}', file: '{args.outfile}')...")
        
        # Connect and import
        with duckdb.connect(str(args.outfile)) as con:
            # We use read_csv_auto but provide some hints for better performance and reliability
            query = f"""
                CREATE TABLE {args.table} AS 
                SELECT * FROM read_csv_auto('{temp_path}', 
                                           header=True, 
                                           delim='\\t', 
                                           quote='', 
                                           escape='', 
                                           sample_size=-1)
            """
            con.execute(query)
            
        print(f"Done! Created '{args.outfile}' with {num_cols} columns.")
        
    except Exception as e:
        print(f"Error occurred: {e}")
        if args.outfile.exists():
            print(f"Cleaning up failed output file '{args.outfile}'...")
            args.outfile.unlink()
        sys.exit(1)
        
    finally:
        if not args.keep_temp and temp_path.exists():
            print(f"Removing temporary file '{temp_path}'...")
            temp_path.unlink()

if __name__ == "__main__":
    main()
