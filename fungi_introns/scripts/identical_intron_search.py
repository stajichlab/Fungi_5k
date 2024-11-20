#!/usr/bin/env python

import pandas as pd
import duckdb
import os

import argparse

def main():
    parser = argparse.ArgumentParser(description='Search for identical introns in fungi')
    parser.add_argument('--db', type=str, default="intronblast.db", help='Path to the database file')
    parser.add_argument('--input', type=str, help='Path to the database CSV/TSV file')
    parser.add_argument('--output', type=str, help='Path to the output file')
    args = parser.parse_args()
    
    if not os.path.exists(args.db):
            with duckdb.connect(args.db) as con:
                df = con.read_csv(args.input, header = False, sep = "\t", 
                                compression="gzip", names = ['query', 'qlength', 'subject', 'slength', 'identity', 'alnlen', 'mismatches', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
                con.register("df_view", df)
                print(df)
                con.execute("CREATE TABLE intron_matches AS SELECT * FROM df_view")        
                con.execute("CREATE INDEX idx_query ON intron_matches(query)")
                con.execute("CREATE INDEX idx_subject ON intron_matches(subject)")
                con.execute("CREATE INDEX idx_identity ON intron_matches(identity)")
                
    with duckdb.connect(args.db) as con:
        df2 = con.sql("SELECT COUNT(*) FROM intron_matches WHERE identity > 95")
        print(df2)

    
if __name__ == "__main__":
    main()