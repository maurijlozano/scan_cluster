#!/usr/bin/env python3

import argparse, sys, os, glob, re
from Bio import SeqIO

def parseArgs():
    '''
    Argument parsing is done.
    '''
    parser = argparse.ArgumentParser(description='Rename tree generated with scan_cluster.')
    inputArgs = parser.add_argument_group('Required files/folders')#Genome with the query cluster, gb file
    inputArgs.add_argument("-g", "--genomeFolder",help="Folder with the genomes in genbank format.", dest="gfolder", action='store', required=True)
    args = parser.parse_args()
    return args
    #start and end cluster coordinates


if __name__ == "__main__":
    args = parseArgs()
    genomeFolder= args.gfolder
    target_genome_files = glob.glob(f'{genomeFolder}/*.gb')
    with open('genomes.csv','w') as tf:
        tf.write('Assembly Accession number,Organism,Taxonomy,Description,Comments\n')
        for gnome in target_genome_files:
            acc = os.path.splitext(os.path.basename(gnome))[0]
            r = list(SeqIO.parse(gnome,'gb'))[0]
            if len(r.annotations['organism']) > 0:
                org = r.annotations['organism']
            else:
                org = ''
            if len(r.annotations['taxonomy']) > 0:
                tax = '; '.join(r.annotations['taxonomy'])
            else:
                tax = ''
            desc = re.sub(',',';',r.description)
            desc = re.sub(' ?[Cc][Oo][Nn][Tt][Ii][Gg].*','',desc)
            desc = re.sub(' ?[Nn][Oo][Dd][Ee].*','',desc)
            desc = re.sub(' ?[Ss][cC][aA][Ff].*','',desc)
            desc = re.sub(' ?[Pp][Ll][aA][Ss].*','',desc)
            desc = re.sub(' ?[Cc][Hh][Rr][Oo].*','',desc)
            if 'comment' in r.annotations.keys():
                if len(r.annotations['comment']) > 0:
                    comments = re.sub('[,]','; ',r.annotations['comment'])
                    comments = re.sub('\n',' ',comments)
            else:
                comments = ''
            tf.write(f'{acc},{org},{tax},{desc},{comments}\n')
                
    print('genomes.csv table was successfully generated..')
    