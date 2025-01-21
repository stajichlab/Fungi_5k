#!/usr/bin/env python

import pandas as pd
import duckdb
import os

import argparse

def main():
    parser = argparse.ArgumentParser(description='Search for identical introns in fungi')
    parser.add_argument('--db', type=str, default="intronblast.db", help='Path to the database file')
    parser.add_argument('--samples', default="../samples.csv", type=str, help='Path to the samples CSV file')
    parser.add_argument('--path', default="blastn_results/*.gz", type=str, help='Path to the samples CSV file')
    parser.add_argument('--input', type=str, help='Path to the database CSV/TSV file')
    parser.add_argument('--output', type=str, help='Path to the output file')
    args = parser.parse_args()
    
    if not os.path.exists(args.db):
            with duckdb.connect(args.db) as con:                
                sf = con.read_csv(args.samples)
                con.register("sf_view", sf)
                con.execute("CREATE TABLE species AS SELECT * FROM sf_view")
                df = con.read_csv(args.input, header = False, sep = "\t", 
                                compression="gzip", names = ['query', 'qlength', 'subject', 'slength', 'identity', 'alnlen', 'mismatches', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
                con.register("df_view", df)
                print(df)
                con.execute("""
                            CREATE TABLE intron_matches AS 
                            SELECT *, substring(query,1,8) as qlocus, substring(subject,1,8) as slocus 
                            FROM df_view WHERE identity >= 90 AND query != subject
                            """)
                
    with duckdb.connect(args.db) as con:
        df2 = con.sql("""
                    SELECT query_species.PHYLUM, query_species.SPECIES, i.qlocus, 
                            subject_species.PHYLUM, subject_species.SPECIES, i.slocus,
                    i.identity,  i.evalue, i.alnlen, qlength, slength,
                    round((qlength/slength),6) as len_ratio,
                    round((alnlen / ( (qlength+slength) / 2)),6) as aln_ratio
                    FROM intron_matches i, species query_species, species subject_species
                    WHERE query_species.LOCUSTAG = i.qlocus AND subject_species.LOCUSTAG = i.slocus AND 
                            identity > 95 and query != subject AND (alnlen / ( (qlength+slength) / 2)) > 0.8 AND
                            query_species.GENUS != subject_species.GENUS ORDER BY identity DESC, alnlen DESC, query_species.SPECIES, qlocus
                    """) # find 95% identical introns
        print(df2)

    
if __name__ == "__main__":
    main()