#!/usr/bin/env python3
"""
Script to dump proteins by Pfam domain from DuckDB functional database.

This script queries the functional database to retrieve protein sequences
that match a specified Pfam domain ID or accession, and outputs them in FASTA format.
"""

import argparse
import duckdb
import sys
import os


def connect_to_database(db_path):
    """Connect to the DuckDB database."""
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database file not found: {db_path}")
    
    return duckdb.connect(db_path,read_only=True)


def build_query(pfam_field, pfam_value):
    """Build the SQL query based on the search field and value."""
    base_query = """
    SELECT 
        species.*,
        gene_proteins.peptide,
        pfam.*
    FROM pfam, gene_proteins, species 
    WHERE {pfam_condition}
        AND species.LOCUSTAG = pfam.species_prefix
        AND pfam.protein_id = gene_proteins.transcript_id
    """.format(pfam_condition=f"{pfam_field} = '{pfam_value}'")
    
    return base_query


def format_fasta_header(row):
    """Format the FASTA header according to specifications."""
    # Extract fields from the row
    species_name = getattr(row, 'SPECIESIN', getattr(row, 'SPECIES', 'Unknown'))
    transcript_id = getattr(row, 'transcript_id', 'Unknown')
    phylum = getattr(row, 'PHYLUM', 'Unknown')
    subphylum = getattr(row, 'SUBPHYLUM', 'Unknown')
    class_name = getattr(row, 'CLASS', 'Unknown')
    pfam_id = getattr(row, 'pfam_id', 'Unknown')
    evalue = getattr(row, 'full_seq_e_value', 'Unknown')
    
    # Calculate sequence length
    peptide = getattr(row, 'peptide', '')
    length = len(peptide) if peptide else 0
    
    header = f">{species_name}_{transcript_id} PHYLUM={phylum} SUBPHYLUM={subphylum} CLASS={class_name} PFAM_ID={pfam_id} PFAM_EVALUE={evalue} LENGTH={length}"
    
    return header


def write_fasta_sequences(cursor, output_file=None):
    """Write sequences in FASTA format to file or stdout."""
    output = open(output_file, 'w') if output_file else sys.stdout
    
    try:
        count = 0
        for row in cursor.fetchall():
            # Convert row to object-like access
            class Row:
                def __init__(self, row_data, column_names):
                    for i, col in enumerate(column_names):
                        setattr(self, col.lower(), row_data[i])
            
            # Get column names from cursor description
            column_names = [desc[0] for desc in cursor.description]
            row_obj = Row(row, column_names)
            
            # Format header and sequence
            header = format_fasta_header(row_obj)
            sequence = getattr(row_obj, 'peptide', '')
            
            if sequence:
                output.write(f"{header}\n")
                # Write sequence in lines of 80 characters
                for i in range(0, len(sequence), 80):
                    output.write(f"{sequence[i:i+80]}\n")
                count += 1
        
        print(f"# Wrote {count} sequences", file=sys.stderr)
        
    finally:
        if output_file:
            output.close()


def main():
    parser = argparse.ArgumentParser(
        description='Dump proteins by Pfam domain from functional database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Search by Pfam ID (default)
  python dump_proteins_by_domain.py --pfam_id Ice_binding
  
  # Search by Pfam accession
  python dump_proteins_by_domain.py --pfam_acc PF11999
  
  # Output to file
  python dump_proteins_by_domain.py --pfam_id Ice_binding -o ice_binding_proteins.fasta
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
        
        # Build and execute query
        query = build_query(pfam_field, pfam_value)
        print(f"# Executing query...", file=sys.stderr)
        cursor = conn.execute(query)
        
        # Write results
        write_fasta_sequences(cursor, args.output)
        
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
