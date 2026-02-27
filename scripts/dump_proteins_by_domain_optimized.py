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


def build_optimized_query(pfam_field, pfam_value, extract_domains=False):
    """Build an optimized SQL query with explicit JOINs."""
    # Use explicit INNER JOINs for better query optimization
    # Select only necessary columns to reduce data transfer
    if extract_domains:
        # Include env_from, env_to, and domain_num for domain extraction
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
            p.full_seq_e_value as evalue,
            p.env_from,
            p.env_to,
            p.domain_num
        FROM pfam p
        INNER JOIN gene_proteins gp ON p.protein_id = gp.transcript_id
        INNER JOIN species s ON s.LOCUSTAG = p.species_prefix
        WHERE p.{pfam_field} LIKE '{pfam_value}%'
        ORDER BY species_name, gp.transcript_id, p.domain_num
        """.format(pfam_field=pfam_field, pfam_value=pfam_value)
    else:
        # Original query without domain extraction fields
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


def format_fasta_header_fast(row_data, col_indices, domain_suffix="", domain_sequence=None, env_from=None, env_to=None):
    """Fast FASTA header formatting using precomputed column indices."""
    species_name = row_data[col_indices['species_name']] or 'Unknown'
    species_name = re.sub(r' ','_',species_name)
    transcript_id = row_data[col_indices['transcript_id']] or 'Unknown'
    phylum = row_data[col_indices['phylum']] or 'Unknown'
    subphylum = row_data[col_indices['subphylum']] or 'Unknown'
    class_name = row_data[col_indices['class']] or 'Unknown'
    pfam_id = row_data[col_indices['pfam_id']] or 'Unknown'
    evalue = row_data[col_indices['evalue']] or 'Unknown'
    
    # Add domain suffix to transcript_id if provided
    transcript_id_with_suffix = f"{transcript_id}{domain_suffix}"
    
    # Use domain sequence length if provided, otherwise use full peptide length
    if domain_sequence is not None:
        length = len(domain_sequence)
    else:
        peptide = row_data[col_indices['peptide']] or ''
        length = len(peptide)
    
    # Build base header
    header = f">{species_name}_{transcript_id_with_suffix} PHYLUM={phylum} SUBPHYLUM={subphylum} CLASS={class_name} PFAM_ID={pfam_id} PFAM_EVALUE={evalue} LENGTH={length}"
    
    # Add domain coordinates if provided
    if env_from is not None and env_to is not None:
        header += f" ENV_FROM={env_from} ENV_TO={env_to}"
    
    return header


def extract_domain_sequence(full_sequence, env_from, env_to):
    """Extract domain sequence using env_from and env_to coordinates."""
    if not full_sequence or env_from is None or env_to is None:
        return ""
    
    # Convert to 0-based indexing (env_from and env_to are 1-based)
    start_pos = max(0, int(env_from) - 1)
    end_pos = min(len(full_sequence), int(env_to))
    
    # Extract the domain sequence
    domain_sequence = full_sequence[start_pos:end_pos]
    return domain_sequence


def write_fasta_batch(cursor, output_file=None, batch_size=1000, extract_domains=False):
    """Write FASTA sequences using batch processing for better I/O performance."""
    output = open(output_file, 'w') if output_file else sys.stdout
    
    try:
        # Get column names and create index mapping for fast access
        column_names = [desc[0].lower() for desc in cursor.description]
        col_indices = {name: idx for idx, name in enumerate(column_names)}
        
        # Verify required columns exist
        required_cols = ['species_name', 'transcript_id', 'phylum', 'subphylum', 
                        'class', 'pfam_id', 'evalue', 'peptide']
        if extract_domains:
            required_cols.extend(['env_from', 'env_to', 'domain_num'])
            
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
        
        if extract_domains:
            # Group rows by transcript_id to handle multiple domains per protein
            from collections import defaultdict
            protein_domains = defaultdict(list)
            
            for row in all_rows:
                transcript_id = row[col_indices['transcript_id']] if col_indices['transcript_id'] < len(row) else 'Unknown'
                protein_domains[transcript_id].append(row)
            
            # Process each protein and its domains
            for i, (transcript_id, domains) in enumerate(protein_domains.items()):
                try:
                    # Count domains for this protein to create proper suffixes
                    domain_count = {}
                    for domain in domains:
                        pfam_id = domain[col_indices['pfam_id']] if col_indices['pfam_id'] < len(domain) else 'Unknown'
                        if pfam_id not in domain_count:
                            domain_count[pfam_id] = 0
                        domain_count[pfam_id] += 1
                        
                        full_sequence = domain[col_indices['peptide']] if col_indices['peptide'] < len(domain) else ''
                        env_from = domain[col_indices['env_from']] if col_indices['env_from'] < len(domain) else None
                        env_to = domain[col_indices['env_to']] if col_indices['env_to'] < len(domain) else None
                        
                        # Extract domain sequence
                        domain_sequence = extract_domain_sequence(full_sequence, env_from, env_to)
                        
                        if domain_sequence:
                            # Create domain suffix
                            domain_suffix = f".domain_{domain_count[pfam_id]}"
                            
                            # Format header with domain suffix, correct length, and coordinates
                            header = format_fasta_header_fast(domain, col_indices, domain_suffix, domain_sequence, env_from, env_to)
                            
                            output_buffer.append(f"{header}\n")
                            # Write sequence in lines of 80 characters
                            for j in range(0, len(domain_sequence), 80):
                                output_buffer.append(f"{domain_sequence[j:j+80]}\n")
                            count += 1
                    
                    # Write in batches for better I/O performance
                    if len(output_buffer) >= batch_size or i == len(protein_domains) - 1:
                        output.write(''.join(output_buffer))
                        output_buffer.clear()
                        
                        # Progress indicator for large datasets
                        if i > 0 and i % 1000 == 0:
                            print(f"# Processed {i+1}/{len(protein_domains)} proteins...", file=sys.stderr)
                            
                except (IndexError, KeyError) as e:
                    print(f"Warning: Skipping protein {transcript_id} due to data issue: {e}", file=sys.stderr)
                    continue
        else:
            # Original behavior for non-domain extraction mode
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


def write_fasta_streaming(cursor, output_file=None, extract_domains=False):
    """Streaming version for very large datasets (lower memory usage)."""
    output = open(output_file, 'w') if output_file else sys.stdout
    
    try:
        # Get column names and create index mapping
        column_names = [desc[0].lower() for desc in cursor.description]
        col_indices = {name: idx for idx, name in enumerate(column_names)}
        
        count = 0
        row_num = 0
        
        if extract_domains:
            print("# Warning: Streaming mode with domain extraction may not group domains properly", file=sys.stderr)
            print("# Consider using batch mode for domain extraction", file=sys.stderr)
            
            # For streaming, we'll process domains individually without grouping
            # This means domain numbering may not be sequential for the same protein
            domain_counters = {}  # Track domain counts per protein
            
            for row in cursor:
                row_num += 1
                try:
                    transcript_id = row[col_indices.get('transcript_id', 0)] or 'Unknown'
                    pfam_id = row[col_indices.get('pfam_id', 0)] or 'Unknown'
                    
                    # Create a key for tracking domain instances
                    protein_domain_key = f"{transcript_id}_{pfam_id}"
                    if protein_domain_key not in domain_counters:
                        domain_counters[protein_domain_key] = 0
                    domain_counters[protein_domain_key] += 1
                    
                    full_sequence = row[col_indices.get('peptide', 0)] or ''
                    env_from = row[col_indices.get('env_from', 0)]
                    env_to = row[col_indices.get('env_to', 0)]
                    
                    # Extract domain sequence
                    domain_sequence = extract_domain_sequence(full_sequence, env_from, env_to)
                    
                    if domain_sequence:
                        # Create domain suffix
                        domain_suffix = f".domain_{domain_counters[protein_domain_key]}"
                        
                        # Format header with domain suffix, correct length, and coordinates
                        header = format_fasta_header_fast(row, col_indices, domain_suffix, domain_sequence, env_from, env_to)
                        
                        output.write(f"{header}\n")
                        # Write sequence in lines of 80 characters
                        for i in range(0, len(domain_sequence), 80):
                            output.write(f"{domain_sequence[i:i+80]}\n")
                        count += 1
                    
                    # Progress indicator
                    if row_num % 10000 == 0:
                        print(f"# Processed {row_num} rows, wrote {count} sequences...", file=sys.stderr)
                        
                except (IndexError, KeyError) as e:
                    print(f"Warning: Skipping row {row_num} due to data issue: {e}", file=sys.stderr)
                    continue
        else:
            # Original behavior for non-domain extraction mode
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
  
  # Extract individual domain sequences from proteins
  python dump_proteins_by_domain_optimized.py --pfam_id Ice_binding --extract_domains -o domains.fasta
  
  # Extract domains with custom batch size for better performance
  python dump_proteins_by_domain_optimized.py --pfam_id Ice_binding --extract_domains --batch_size 2000
  
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
        '--extract_domains',
        action='store_true',
        help='Extract individual domain sequences using env_from/env_to coordinates. Each domain instance will be output as a separate sequence with suffix .domain_N (where N is the domain number for that protein)'
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
    
    if args.extract_domains:
        print(f"# Domain extraction mode enabled", file=sys.stderr)
        print(f"# Each domain instance will be output as a separate sequence", file=sys.stderr)
        if args.streaming:
            print(f"# WARNING: Streaming mode with domain extraction is experimental", file=sys.stderr)
            print(f"# Consider using batch mode for better reliability", file=sys.stderr)
    
    # Initialize connection variable
    conn = None
    
    try:
        # Connect to database
        print(f"# Connecting to database: {args.database}", file=sys.stderr)
        conn = connect_to_database(args.database)
        
        # Build and execute optimized query
        query = build_optimized_query(pfam_field, pfam_value, args.extract_domains)
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
            write_fasta_streaming(cursor, args.output, args.extract_domains)
        else:
            write_fasta_batch(cursor, args.output, args.batch_size, args.extract_domains)
        
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
