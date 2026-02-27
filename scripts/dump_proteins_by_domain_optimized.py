#!/usr/bin/env python3
"""
Optimized version of dump_proteins_by_domain.py script.

This script includes several performance optimizations:
1. Improved SQL query with explicit JOINs for better optimization
2. Batch processing for I/O operations
3. More efficient string handling and memory usage
4. Optional streaming mode for very large datasets
"""

import argparse
import duckdb
import sys
import os
import re


def connect_to_database(db_path,threads=4,mem=2):
    """Connect to the DuckDB database with optimizations."""
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database file not found: {db_path}")
    
    conn = duckdb.connect(db_path,read_only=True)
    
    # Set DuckDB optimizations
    conn.execute(f"SET threads={threads}")  # Use multiple threads
    conn.execute(f"SET memory_limit='{mem}GB'")  # Allow more memory usage
    
    return conn


def build_optimized_query(pfam_field, pfam_value):
    """Build an optimized SQL query with explicit JOINs."""
    # Use explicit INNER JOINs for better query optimization
    # Select only necessary columns to reduce data transfer
    optimized_query = """
    SELECT 
        s.SPECIESIN as species_name,
        s.phylum,
        s.subphylum,
        s.CLASS as class,
        gp.peptide,
        gp.transcript_id,
        p.pfam_id,
        p.pfam_acc,
        p.full_seq_e_value as evalue
    FROM pfam p
    INNER JOIN gene_proteins gp ON p.protein_id = gp.transcript_id
    INNER JOIN species s ON s.LOCUSTAG = p.species_prefix
    WHERE p.{pfam_field} LIKE '{pfam_value}%'
    ORDER BY species_name, gp.transcript_id
    """.format(pfam_field=pfam_field, pfam_value=pfam_value)
    
    return optimized_query


def format_fasta_header_fast(row_data, col_indices):
    """Fast FASTA header formatting using precomputed column indices."""
    species_name = row_data[col_indices['species_name']] or 'Unknown'
    species_name = re.sub(r' ','_',species_name)
    transcript_id = row_data[col_indices['transcript_id']] or 'Unknown'
    phylum = row_data[col_indices['phylum']] or 'Unknown'
    subphylum = row_data[col_indices['subphylum']] or 'Unknown'
    class_name = row_data[col_indices['class']] or 'Unknown'
    pfam_id = row_data[col_indices['pfam_id']] or 'Unknown'
    evalue = row_data[col_indices['evalue']] or 'Unknown'
    peptide = row_data[col_indices['peptide']] or ''
    
    length = len(peptide) if peptide else 0
    
    return f">{species_name}_{transcript_id} PHYLUM={phylum} SUBPHYLUM={subphylum} CLASS={class_name} PFAM_ID={pfam_id} PFAM_EVALUE={evalue} LENGTH={length}"


def write_fasta_batch(cursor, output_file=None, batch_size=1000):
    """Write FASTA sequences using batch processing for better I/O performance."""
    output = open(output_file, 'w') if output_file else sys.stdout
    
    try:
        # Get column names and create index mapping for fast access
        column_names = [desc[0].lower() for desc in cursor.description]
        col_indices = {name: idx for idx, name in enumerate(column_names)}
        
        # Verify required columns exist
        required_cols = ['species_name', 'transcript_id', 'phylum', 'subphylum', 
                        'class', 'pfam_id', 'evalue', 'peptide']
        missing_cols = [col for col in required_cols if col not in col_indices]
        if missing_cols:
            print(f"Warning: Missing columns in query result: {missing_cols}", file=sys.stderr)
            # Use safe defaults for missing columns
            for col in missing_cols:
                col_indices[col] = 0  # Will default to first column or cause IndexError
        
        count = 0
        output_buffer = []
        
        # Process all rows
        all_rows = cursor.fetchall()
        print(f"# Processing {len(all_rows)} database rows", file=sys.stderr)
        
        for i, row in enumerate(all_rows):
            try:
                # Format header and sequence
                header = format_fasta_header_fast(row, col_indices)
                sequence = row[col_indices['peptide']] if col_indices['peptide'] < len(row) else ''
                
                if sequence:
                    output_buffer.append(f"{header}\n")
                    # Write sequence in lines of 80 characters
                    for j in range(0, len(sequence), 80):
                        output_buffer.append(f"{sequence[j:j+80]}\n")
                    count += 1
                
                # Write in batches for better I/O performance
                if len(output_buffer) >= batch_size or i == len(all_rows) - 1:
                    output.write(''.join(output_buffer))
                    output_buffer.clear()
                    
                    # Progress indicator for large datasets
                    if i > 0 and i % 10000 == 0:
                        print(f"# Processed {i+1}/{len(all_rows)} rows...", file=sys.stderr)
                        
            except (IndexError, KeyError) as e:
                print(f"Warning: Skipping row {i+1} due to data issue: {e}", file=sys.stderr)
                continue
        
        print(f"# Wrote {count} sequences", file=sys.stderr)
        
    finally:
        if output_file:
            output.close()


def write_fasta_streaming(cursor, output_file=None):
    """Streaming version for very large datasets (lower memory usage)."""
    output = open(output_file, 'w') if output_file else sys.stdout
    
    try:
        # Get column names and create index mapping
        column_names = [desc[0].lower() for desc in cursor.description]
        col_indices = {name: idx for idx, name in enumerate(column_names)}
        
        count = 0
        row_num = 0
        
        # Process one row at a time (memory efficient)
        for row in cursor:
            row_num += 1
            try:
                # Format header and sequence
                header = format_fasta_header_fast(row, col_indices)
                sequence = row[col_indices.get('peptide', 0)] or ''
                
                if sequence:
                    output.write(f"{header}\n")
                    # Write sequence in lines of 80 characters
                    for i in range(0, len(sequence), 80):
                        output.write(f"{sequence[i:i+80]}\n")
                    count += 1
                
                # Progress indicator
                if row_num % 10000 == 0:
                    print(f"# Processed {row_num} rows, wrote {count} sequences...", file=sys.stderr)
                    
            except (IndexError, KeyError) as e:
                print(f"Warning: Skipping row {row_num} due to data issue: {e}", file=sys.stderr)
                continue
        
        print(f"# Wrote {count} sequences (streaming mode)", file=sys.stderr)
        
    finally:
        if output_file:
            output.close()


def main():
    parser = argparse.ArgumentParser(
        description='Dump proteins by Pfam domain from functional database (optimized version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Search by Pfam ID (default)
  python dump_proteins_by_domain_optimized.py --pfam_id Ice_binding
  
  # Search by Pfam accession
  python dump_proteins_by_domain_optimized.py --pfam_acc PF11999
  
  # Output to file with streaming (for very large results)
  python dump_proteins_by_domain_optimized.py --pfam_id Ice_binding -o results.fasta --streaming
  
  # Use larger batch size for better performance
  python dump_proteins_by_domain_optimized.py --pfam_id Ice_binding --batch_size 5000
        """
    )
    
    parser.add_argument(
        '--pfam_id',
        default='Ice_binding',
        help='Pfam domain ID to search for (default: Ice_binding)'
    )
    
    parser.add_argument(
        '--pfam_acc',
        help='Pfam accession to search for (e.g., PF11999)'
    )
    
    parser.add_argument(
        '--database',
        default='functionalDB/function.duckdb',
        help='Path to DuckDB database file (default: functionalDB/function.duckdb)'
    )
    
    parser.add_argument(
        '-o', '--output',
        help='Output file (default: stdout)'
    )
    
    parser.add_argument(
        '--streaming',
        action='store_true',
        help='Use streaming mode (lower memory usage, good for very large datasets)'
    )
    
    parser.add_argument(
        '--batch_size',
        type=int,
        default=1000,
        help='Batch size for I/O operations (default: 1000)'
    )
    
    args = parser.parse_args()
    
    # Determine which field to search on
    if args.pfam_acc:
        pfam_field = 'pfam_acc'
        pfam_value = args.pfam_acc
        print(f"# Searching for Pfam accession: {pfam_value}", file=sys.stderr)
    else:
        pfam_field = 'pfam_id'
        pfam_value = args.pfam_id
        print(f"# Searching for Pfam ID: {pfam_value}", file=sys.stderr)
    
    # Initialize connection variable
    conn = None
    
    try:
        # Connect to database
        print(f"# Connecting to database: {args.database}", file=sys.stderr)
        conn = connect_to_database(args.database)
        
        # Build and execute optimized query
        query = build_optimized_query(pfam_field, pfam_value)
        print(f"# Executing optimized query...", file=sys.stderr)
        
        # Show query plan for debugging (optional)
        if os.getenv('DEBUG'):
            print("# Query plan:", file=sys.stderr)
            explain_result = conn.execute(f"EXPLAIN {query}").fetchall()
            for row in explain_result:
                print(f"#   {row[0]}", file=sys.stderr)
        
        cursor = conn.execute(query)
        
        # Write results using selected method
        if args.streaming:
            write_fasta_streaming(cursor, args.output)
        else:
            write_fasta_batch(cursor, args.output, args.batch_size)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    finally:
        if conn is not None:
            try:
                conn.close()
            except:
                pass


if __name__ == '__main__':
    main()
