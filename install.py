#!/usr/bin/env python3
print('Installing the requirements for Scan cluster...')
import os,sys
from sys import platform
from shutil import which
import subprocess

os.system('pip install -r requirements.txt')

if not platform.startswith("linux"):
    print('This program only works in linux systems...\nFor windows or Mac OS you will have to manually install Blast and HMMER.')

def is_insatalled(program):
    """Check whether `program` is on PATH and marked as executable."""
    return (which(program) != None)

if not is_insatalled('blastp'):
    install_blast = input('Do you want to install NCBI Blast? [Yes/No]')
    if (install_blast.upper == 'YES') or (install_blast.upper == 'Y'): 
        os.system('sudo apt install ncbi-blast+')
else:
    print(f'Blast is installed: {which("blastp")}')

if not is_insatalled('hmmsearch'):
    install_hmmer = input('Do you want to install HMMER? [Yes/No]')
    if (install_hmmer.upper == 'YES') or (install_hmmer.upper == 'Y'): 
        os.system('sudo apt install hmmer')
else:
    print(f'HMMER is installed: {which("hmmsearch")}')

if not is_insatalled('mafft'):
    install_mafft = input('Do you want to install mafft? [Yes/No]')
    if (install_mafft.upper == 'YES') or (install_mafft.upper == 'Y'): 
        os.system('sudo apt install mafft')
else:
    print(f'mafft is installed: {which("mafft")}')

if os.path.exists('scan_cluster.py'):
    os.system('chmod +x scan_cluster.py')
    print(f'\n\nTesting scan_cluster')
    os.system('./scan_cluster.py -h')
else:
    print('scan_cluster.py not found in current folder.')
