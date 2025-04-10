#!/usr/bin/env python3
import pandas as pd
import sys, os

args = sys.argv
if len(args) != 4:
    sys.exit('Error, tsv table file with headers, accession number field name, and output folder required...')

# args =['','ncbi_dataset.tsv']
table = pd.read_csv(args[1], sep='\t')
#accession_column_name = 'Assembly Release Date'
accession_column_name = args[2]
if not accession_column_name in table.columns:
    sys.exit('Error: Accession number field not found in table...')

output_folder = args[3]
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

table = table.sort_values([accession_column_name], ascending=False)
acc = [an for an in table[accession_column_name] ]
for rs in acc:
    ofile = os.path.join(output_folder,rs+".gb")
    if (not os.path.exists(ofile)) or (os.path.getsize(ofile) < 1000 ):
        os.system(f'datasets download genome accession {rs} --include gbff --filename {rs}.zip')
        os.system(f'unzip -p {rs}.zip ncbi_dataset/data/{rs}/genomic.gbff > {ofile}')
        os.remove(f'{rs}.zip')
    else:
        print(f'{rs} already downloaded...')
