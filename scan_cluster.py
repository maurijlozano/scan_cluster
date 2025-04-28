#!/usr/bin/env python3
#requirements Blast >2.13
NAME = "Scan cluster"
VERSION = "1.0"
REF = "\n   Not published"
GITHUB="https://github.com/maurijlozano"
BLAST_DBS = ['nr', 'refseq_select_prot', 'refseq_protein', 'SMARTBLAST/landmark', 'swissprot', 'pataa', 'pdb', 'env_nr', 'tsa_nr', 'nr_cluster_seq']

#Imports
#Python modules
import argparse, sys, os, subprocess, re, glob
from datetime import datetime
from copy import deepcopy
from datetime import datetime

#Third party modules
from Bio import SeqIO
from Bio.Phylo import TreeConstruction
from Bio import Phylo
import pandas as pd

#Functions
def parseArgs():
	'''
	Argument parsing is done.
	'''
	parser = argparse.ArgumentParser(description='This program was designed to identify orthologous genes clusters.')
	inputArgs = parser.add_argument_group('Query cluster')#Genome with the query cluster, gb file
	inputArgs.add_argument("-Q", "--QueryFile",help="Query genome in Genbank format.", dest="qfile", action='store', required=False)
	inputArgs.add_argument("-R", "--Replicon_ID",help="Required for draft and multireplicon genomes.", dest="repid", action='store', required=False)
	#start and end cluster coordinates
	inputArgs.add_argument("-s", "--cluster_start",help="Genome cluster start location.", dest="cstart", action='store', required=False)
	inputArgs.add_argument("-e", "--cluster_end",help="Genome cluster end location.", dest="cend", action='store', required=False)
	#add cluster_gb_file --> gb file containing only the cluster region
	gbCluster = parser.add_argument_group('To search or a cluster provided in Genbank format')
	gbCluster.add_argument("-q", "--QueryCluster",help="Query cluster in Genbank format.", dest="qcluster", action='store', required=False)
	gbCluster.add_argument("--Reference",help="Specifies the genome Gb file to be used as reference.", dest="refgenome", action='store', required=False)
	hmmCluster = parser.add_argument_group('To search or a cluster with a defined set of protein HMMs')
	hmmCluster.add_argument("-f", "--hmm_folder",help="Folder containing the HMM profiles for the proteins to include in the cluster.", dest="hmm_folder", action='store')
	#Subject genome
	targetinput = parser.add_argument_group('Target Genomes')
	targetinput.add_argument("-S", "--SubjectFile",help="Subject genome in Genbank format.", dest="sfile", action='store', required=False)
	targetinput.add_argument("-F", "--SubjectFolder",help="Folder containing the subject genomes in Genbank format.", dest="sfolder", action='store', required=False)
	targetinput.add_argument("-E", "--Genbank file extension",help="Genbank file extension. Default   .gb  ", dest="gbext", action='store', required=False, default='gb')
	#hmm folder
	outopt = parser.add_argument_group('Output Folder')
	outopt.add_argument("-o", "--Results_folder",help="Results folder name.", dest="res_folder", action='store', default = 'Results')
	outopt.add_argument("--overwrite", help="Overwrite previous results.", dest="overwrite", action='store_true')
	#options
	runmode = parser.add_argument_group('Running mode')
	runmode.add_argument("--only_blastp",help="Runs using only blast for the identification of homolog proteins.", dest="only_blastp", action='store_true', required=False)
	#add maximum number of genes not belonging to cluster to end cluster definition
	clusterOPts = parser.add_argument_group('Cluster definition arguments')
	clusterOPts.add_argument("-n", "--n_prots_between",help="Maximum number of proteins allowed between two consecutive genes in the query cluster. Default = half of proteins in the cluster", dest="prots_between", action='store')
	clusterOPts.add_argument("-M", "--max_alien_prots",help="Maximum number of proteins in the target cluster that are not present in the query cluster. Default = not limited (Number of proteins in cluster * 3).", dest="max_alien_prots", action='store')
	clusterOPts.add_argument('-m', "--min_target_prots",help="Minimum of query proteins required to be found in target cluster. Default=3.)", dest="min_target_prots", action='store', default=3)
	clusterOPts.add_argument("--min_cluster_coverage",help="Minimum of cluster coverage, proportion. Default=.5. The program will use as minimum half of the query proteins. If you are running only with HMMs, this value should be the fraction of the HMM required in a cluster.)", dest="min_cluster_coverage", action='store', default=0.5)
	clusterOPts.add_argument('-g', "--gap_penalty",help="Gap penalty for cluster alignment. Default = 10", dest="gap", action='store', default=10)
	clusterOPts.add_argument("--mismatch_score",help="Mismatch score for cluster alignment. Alignment of genes that are not orthologs are penalized. Default = 20", dest="mismatch", action='store', default=20)
	#Optional arguments for blast and HMM
	searchOpt = parser.add_argument_group('Blast and HMMSearch options')
	searchOpt.add_argument("--local_blast_db",help="A local blastp database generated with makeblastdb program... <Folder name>", dest="local_blast_db", action='store', default='')
	searchOpt.add_argument("--Generate_local_db",help="A local blastp database will be generated from the proteome of all the analyzed subject sequences...", dest="local_blast_db_subject", action='store_true')
	searchOpt.add_argument("--Blast_DB",help="Database for remote blastp, used to retrieve homologs for HMM generation. Default = refseq_protein. Available: nr, refseq_select, refseq_protein, landmark, swissprot, pataa, pdb, env_nr, tsa_nr", dest="blastp_database", action='store', default='refseq_protein')
	#
	searchOpt.add_argument("--Blast_evalue",help="E-value cut-off for remote blastp, used to retrieve homologs for HMM generation.", dest="evalue", action='store', default=0.00001)
	searchOpt.add_argument("--Blast_max_targets",help="Maxímum number of targets for Blastp search. Default=250", dest="max_target", action='store', default=250)
	searchOpt.add_argument("--Blast_qcov",help="Query coverage percent for Blastp search. Default=45", dest="qcov", action='store', default=45)
	searchOpt.add_argument("--Blast_scov",help="Subject coverage percent for Blastp search. Default=45", dest="scov", action='store', default=45)
	searchOpt.add_argument("--hmm_evalue",help="E-value cut-off for hmmsearch. Default = 0.00001", dest="hmm_evalue", action='store', default=0.00001)
	searchOpt.add_argument("--hmm_cover",help="HMM coverage cut-off for hmmsearch. Default = 45", dest="hmm_cover", action='store', default=45)	
	searchOpt.add_argument("--ali_cover",help="Alignment coverage cut-off for hmmsearch. Default = 45", dest="ali_cover", action='store', default=45)
	#
	args = parser.parse_args()
	return args

def printSoftName():
	print("\n\n\n")
	print("   ************************************************")
	print("   *****   Scan cluster ")
	print("   *****   Version: "+str(VERSION))
	print("   *****   Developed by Mauricio J. Lozano")
	print("   *****   github.com/maurijlozano")
	print("   ************************************************")
	print("   Please cite: "+REF)
	print("   Downloaded from: "+GITHUB)
	print("\n\n\n")
	now = datetime.now()
	dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
	return(f'Date and time: {dt_string}')

def write_to_log(log_file,message):
	with open(log_file, 'a') as log:
		log.write(f'{message}\n')
	pass

def AnnotatedGenome(genbankFile):
	features = 0
	for record in SeqIO.parse(genbankFile, "genbank"):
		if record.features:
			if len(record.features) <= 5:
				features += 0
			else:
				features += len(record.features)
	if features == 0:
		write_to_log(log_file,f'Error: {genbankFile} is not annotated. Please provide a Genbank file with annotations.')
		print(f'Error: {genbankFile} is not annotated. Please provide a Genbank file with annotations.')
		return(False)
	else:
		return(True)
			
#Genbankfile to dict
def genbankToDict(genbankFile):
	genbankDict = {}
	description = {}
	count=1
	for record in SeqIO.parse(genbankFile, "genbank"):
		featuresDict = {}
		if record.features:
			for feat in record.features:
				if feat.type == "CDS":
					if not 'pseudo' in feat.qualifiers:
						UID = "UID"+str(count)
						proteinSeq = feat.qualifiers['translation'][0]
						if "locus_tag" in feat.qualifiers.keys():
							locusTag = feat.qualifiers["locus_tag"][0]
						elif "gene" in feat.qualifiers.keys():
							locusTag = feat.qualifiers["gene"][0]
						else:
							locusTag = "ND"
						if 'protein_id' in feat.qualifiers:
							proteinID = feat.qualifiers['protein_id'][0]
						else:
							proteinID = "ND"
						featuresDict[UID] = [locusTag,proteinID,proteinSeq,int(feat.location.start), int(feat.location.end), int(feat.location.strand), UID]
						count+=1
				if feat.type == 'source':
					if 'strain' in feat.qualifiers.keys():
						strain = feat.qualifiers['strain'][0]
						organism = feat.qualifiers['organism'][0] + ' (' + strain +')'
					elif 'isolate' in feat.qualifiers.keys():
						strain = feat.qualifiers['isolate'][0]
						organism = feat.qualifiers['organism'][0] + ' (' + strain +')'
					else:
						organism = feat.qualifiers['organism'][0]
		description[record.id] = organism
		genbankDict[record.id] = featuresDict
	return(genbankDict,description)

#Extract translation for all proteins (proteome) from genbank file
def extractProteome(genbankDict,sres_folder,sgenome_name):
	genbankFileBaseName = os.path.join(sres_folder,sgenome_name)
	proteomeFile = genbankFileBaseName+".faa"
	if not os.path.exists(proteomeFile):
		with open(proteomeFile, "w+") as f:
			for key,value in genbankDict.items():
				for k,val in genbankDict[key].items():
					f.write(">"+key+"|"+val[6]+"|"+val[0]+"\n"+str(val[2])+"\n")
	
	return(proteomeFile)

#Search cheA hmm in proteome
def hmmsearchProtein(proteomeFile,hmm_file,sres_folder,hmm_evalue):
	resFile = os.path.join(sres_folder,"hmmsearch.txt")
	subprocess.call(['hmmsearch', "--domtblout", resFile , hmm_file , proteomeFile], stdout=subprocess.DEVNULL)
	tableHeaders = ['target name', 'taccession', 'tlen', 'query name', 'qaccession', 'qlen', 'E-value', 'score', 'bias', '#', 'of', 'c-Evalue', 'i-Evalue', 'dom score', 'dom  bias', 'hmm from', 'hmm to', 'ali from', 'ali to', 'env from', 'env to', 'acc', 'description of target']
	try:
		hmmHits = pd.read_fwf(resFile, comment='#', sep='\\s+', header=None)
	except:
		print(f'---> No homologs found. Empty hmmsearch.txt file...')
		return(pd.DataFrame(columns = tableHeaders))
	if len(hmmHits) > 0:
		renameDict = { old:new for old,new in zip(hmmHits.columns[0:len(tableHeaders)],tableHeaders)}
		hmmHits = hmmHits.rename(columns=renameDict)
		#old vers
		#hmmHits = pd.read_table(resFile, names=tableHeaders, comment='#', sep='\\s+').sort_values(["target name","E-value"])
		hmmHits = hmmHits.sort_values("E-value")
		hmmHitsPerGene = hmmHits.groupby(['target name'], as_index=False, ).first().sort_values("E-value")
		hmmHitsPerGene = hmmHitsPerGene[hmmHitsPerGene['E-value'] < hmm_evalue]
		hmmHitsPerGene['hmm cover'] = [(row['hmm to']-row['hmm from'])/row['qlen'] for _,row in hmmHitsPerGene.iterrows()]
		hmmHitsPerGene['ali cover'] = [(row['ali to']-row['ali from'])/row['tlen'] for _,row in hmmHitsPerGene.iterrows()]
		hmmHitsPerGene = hmmHitsPerGene[hmmHitsPerGene['hmm cover'] >= hmm_cover]
		hmmHitsPerGene = hmmHitsPerGene[hmmHitsPerGene['ali cover'] >= ali_cover]
		hmmHitsPerGene = hmmHitsPerGene.reset_index(drop=True)
		return(hmmHitsPerGene)
	else:
		return(pd.DataFrame(columns = tableHeaders))

#gets the clusters with gene orientation and score
def get_oriented_cluster_with_score(region_file, hit_table,proteins):
	cluster_structure = []
	with open(proteins, 'a') as pfile:
		for r in SeqIO.parse(region_file,'gb'):
			rid = r.id
			for feat in r.features:
				if (feat.type == 'CDS') and ('translation' in feat.qualifiers.keys()):
					size = int(feat.location.end) - int(feat.location.start)
					strand = feat.location.strand
					if "locus_tag" in feat.qualifiers.keys():
						locusTag = feat.qualifiers["locus_tag"][0]
					elif "gene" in feat.qualifiers.keys():
						locusTag = feat.qualifiers["gene"][0]
					else:
						locusTag = "ND"
					repl_locusTag = f'{rid}|{locusTag}'
					hit = hit_table[hit_table['replicon_locus_tag'] == repl_locusTag]
					if len(hit) > 0:
						query = hit['query name'].iloc[0]
					else:
						query = ''
					cluster_structure.append((repl_locusTag,query,size,strand))
					pfile.write(f'>{rid}|{locusTag}\n{feat.qualifiers["translation"][0]}\n')
	return(cluster_structure)

#Extract translation of CDS in cluster region
def extract_cluster(genbankFile,replicon_id,cstart,cend,res_folder):
	query_cluster_faa = os.path.join(res_folder,'query_cluster.faa')
	description = os.path.splitext(os.path.basename(genbankFile))[0]
	with open(query_cluster_faa,'w') as faa:
		id_found = False
		if replicon_id:
			for r in SeqIO.parse(genbankFile,'gb'):
				if r.id == replicon_id:
					query_cluster = r[cstart:cend]
					id_found = True
			if not id_found:
				write_to_log(log_file,f'Error: Replicon/contig ID ({replicon_id}) not found...')
				sys.exit(f'Error: Replicon/contig ID ({replicon_id}) not found...')
		else:
			query_clusters = list(SeqIO.parse(genbankFile,'gb'))
			if len(query_clusters) > 1:
				write_to_log(log_file,'Error: query cluster should contain only one Genbank record...\nThe first record will be used...')
				print('Error: query cluster should contain only one Genbank record...\nThe first record will be used...')
			query_cluster = query_clusters[0]
		if 'features' in vars(query_cluster):
			for feat in query_cluster.features:
				locusTag = ''
				proteinID = ''
				if feat.type == "CDS":
					if not 'pseudo' in feat.qualifiers:
						if not 'translation' in feat.qualifiers.keys():
							continue
						else:
							proteinSeq = feat.qualifiers['translation'][0]
						if "locus_tag" in feat.qualifiers.keys():
							locusTag = feat.qualifiers["locus_tag"][0]
						elif "gene" in feat.qualifiers.keys():
							locusTag = feat.qualifiers["gene"][0]
						if 'protein_id' in feat.qualifiers:
							proteinID = feat.qualifiers['protein_id'][0]
						translation = feat.qualifiers['translation'][0]
						faa.write(f'>{locusTag}|{proteinID}|{description}\n{translation}\n')
		else:
			write_to_log(log_file,f'No CDS were found in the query cluster...')
			sys.exit(f'No CDS were found in the query cluster...')
	return(query_cluster_faa)

#Function that makes a remote blast search against NCBI refseq_proteins, then
# makes a MSA with MAFFT L-INS-I, and generates an HMM
def generate_HMM(query_cluster_faa,evalue,blastp_database,qcov,scov,res_folder,local_blast_db,hmm_folder):
	HMMs = []
	blast_out = os.path.join(res_folder,'blast.res')
	temp_file = os.path.join(res_folder,'temp.faa')
	for protein in SeqIO.parse(query_cluster_faa,'fasta'):
		ltag = protein.id.split('|')[0]
		hmm = os.path.join(hmm_folder,f'{ltag}.hmm')
		if (os.path.exists(hmm)) and (os.path.getsize(hmm) > 0):
			HMMs.append(hmm)
			continue
		SeqIO.write(protein,temp_file,'fasta')
		if local_blast_db == '':	
			blast_command = f"blastp -db {blastp_database} -query {temp_file} -remote -out {blast_out} -max_target_seqs {max_target} -word_size 5 -evalue {evalue} -outfmt '6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore sseq' 1>> {log_file}"
			write_to_log(log_file,f'Running blastp in remote mode to retrieve {ltag} homologs for HMM construction.')
			print(f'> Running blastp in remote mode to retrieve {ltag} homologs for HMM construction.')
		else:
			blast_command = f"blastp -db {local_blast_db} -query {temp_file} -out {blast_out} -max_target_seqs {max_target} -word_size 5 -evalue {evalue} -outfmt '6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore sseq' 1>> {log_file}"
			write_to_log(log_file,f'Running blastp in local mode to retrieve {ltag} homologs for HMM construction.')
			print(f'> Running blastp in local mode to retrieve {ltag} homologs for HMM construction.')
		try:
			os.system(blast_command)
		except:
			message = f'Error: Blast search failed. {blast_command}'
			write_to_log(log_file,message)
			print(message)
		#Check blast results
		if os.path.getsize(blast_out) == 0:
			write_to_log(log_file,'Error: No homologs were found in the local/remote {local_blast_db} database... Please try again or run in blast mode without HMM generation...')
			print(f'Error: No homologs were found for {protein.id} in the local/remote {local_blast_db} database... Please try again or run in blast mode without HMM generation...')
			continue
		#read blast resutls
		blast_fields = ['qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'evalue', 'bitscore', 'sseq']
		blast_tabl = pd.read_table(blast_out, names = blast_fields)
		blast_tabl['qcov'] = [ (row['qend']-row['qstart'])/row['qlen'] for _,row in blast_tabl.iterrows()]
		blast_tabl['scov'] = [ (row['send']-row['sstart'])/row['slen'] for _,row in blast_tabl.iterrows()]
		blast_tabl = blast_tabl[(blast_tabl['qcov'] >= qcov ) & (blast_tabl['scov'] >= scov )]
		#read blast results for each protein
		blast_tabl.sort_values(['pident'])
		if len(blast_tabl) > 100:
			pick_one_every = len(blast_tabl)//100
			indices_to_pick = [ i for i in range(0,len(blast_tabl),pick_one_every) ][0:100]
			subtable_pick100 = blast_tabl.iloc[indices_to_pick,]
		else:
			subtable_pick100 = blast_tabl
		#
		blast_fasta = os.path.join(res_folder,f'blast_{ltag}.fasta')
		blast_fasta_aln = os.path.join(res_folder,f'blast_{ltag}_aln.fasta')
		with open(blast_fasta,'w') as fasta:
			for _,row in subtable_pick100.iterrows():
				fasta.write(f">{row['sseqid']}\n{row['sseq']}\n")
			if len(subtable_pick100) == 1:
				fasta.write(f">{row['sseqid']}\n{row['sseq']}\n")
				write_to_log(log_file,f'Only one hit for {protein.id}...')
				print(f'--> Only one hit for {protein.id}...')
			elif len(subtable_pick100) == 0:
				write_to_log(log_file,f'No hits for {protein.id}...')
				print(f'--> No hits for {protein.id}...')
		#make alignment with mafft
		write_to_log(log_file,f'Making {ltag} MSA with Mafft...')
		print(f'--> Making {ltag} MSA with Mafft...')
		os.system(f'mafft --quiet --localpair --maxiterate 1000 {blast_fasta} > {blast_fasta_aln}')
		#make HMM with hmmbuild
		os.system(f'hmmbuild -n {ltag} {hmm} {blast_fasta_aln} 1>> {log_file}')
		HMMs.append(hmm)
	if os.path.exists(temp_file):
		os.remove(temp_file)
	if os.path.exists(blast_out):
		os.remove(blast_out)
	if len(HMMs) > 0:
		return(HMMs)
	else:
		sys.exit(f'No homologs found, could not create the HMMs. Try running in blast only mode.')

#finds by composition, consecutive genes, but not synteny or orientation of the genes
def find_clusters(protHits,nprot,prots_between,max_alien_prots,min_target_prots,min_cluster_coverage):
	#get contiguous
	protHits.loc[:,'taccession'] = [ n.split('|')[0] for n in protHits['target name']]
	protHits.loc[:,'tuid'] = [int(re.sub('UID','',n.split('|')[1])) for n in protHits['target name']]
	clusters = []
	contigsWithHit = protHits['taccession'].unique()
	for cwh in contigsWithHit:
		protHits_ctg = protHits[protHits['taccession'] == cwh]
		prot_dict = { row['target name'].split('|')[1]:row['query name'] for _,row in protHits_ctg.iterrows()}	
		prot_numbers = [ int(t.split('|')[1][3:]) for t in protHits_ctg['target name']]
		prot_numbers.sort()		
		prots_in_cluster = []
		for i in range(len(prot_numbers)-1):
			if len(prots_in_cluster) == 0:
				prots_in_cluster.append(prot_numbers[i])
			if prot_numbers[i+1]-prot_numbers[i] <= prots_between:
				prots_in_cluster.append(prot_numbers[i+1])
				if i == len(prot_numbers)-2:
					clusters.append(prots_in_cluster)	
			else:
				clusters.append(prots_in_cluster)
				prots_in_cluster = []
	#
	clusters = [ clust for clust in clusters if len(clust) >= min_target_prots]
	#
	validated_clusters = []
	validated_clusters_UIDS = []
	for cluster in clusters:
		nprots_found_in_cluster = len(cluster)
		cluster_coverage = (nprots_found_in_cluster)/nprot
		if cluster_coverage < min_cluster_coverage:
			print(f'--> Cluster {cluster} discarded due to low cluster coverage: {cluster_coverage} < {min_cluster_coverage}')
			continue
		else:
			cluster_structure = list(range(min(cluster),max(cluster)+1))
			cluster_structure_homologs = [ prot_dict[f'UID{number}'] if f'UID{number}' in prot_dict.keys() else f'UID{number}' for number in cluster_structure]
			cluster_structure_UIDS = [ f'UID{number}' for number in cluster_structure]
			total_prots_in_cluster = len(cluster_structure)
			alien_prots = total_prots_in_cluster-nprots_found_in_cluster
			if alien_prots >= max_alien_prots:
				print(f'--> Cluster {cluster} discarded due to an excessive number of alien proteins: {alien_prots} > {max_alien_prots}')
				continue
			validated_clusters.append(cluster_structure_homologs)
			validated_clusters_UIDS.append(cluster_structure_UIDS)
	#
	protHits.loc[:,'cluster'] = ''
	for i,c in enumerate(validated_clusters_UIDS):
		protHits.loc[:,'cluster'] = [ i if row['target name'].split('|')[1] in c else row['cluster'] for _,row in protHits.iterrows() ]
	#
	protHits = protHits[protHits['cluster'] != '']
	return(validated_clusters,validated_clusters_UIDS,protHits)

def scan_cluster_HMMER(HMMs, nprot, prots_between,max_alien_prots, sgenbankFile, res_folder, min_target_prots,min_cluster_coverage,hmm_evalue):
	if not os.path.exists(sgenbankFile):
		write_to_log(log_file,"Error: The genome genbank file could not be found.")
		sys.exit("Error: The genome genbank file could not be found.")
	try:
		sgenome_name = os.path.splitext(os.path.basename(sgenbankFile))[0]
		sres_folder = os.path.join(res_folder,sgenome_name)
		if not os.path.exists(sres_folder):
			os.mkdir(sres_folder)
		sgenbankDict, description = genbankToDict(sgenbankFile)
		proteomeFile = extractProteome(sgenbankDict,sres_folder,sgenome_name)
	except:
		write_to_log(log_file,f"Error: Is {sgenome_name} genome in genbank format?")
		sys.exit(f"Error: Is {sgenome_name} genome in genbank format?")
	#search for HMMs
	hmm_file = 'temp.hmm'
	os.system(f'cat {" ".join(HMMs)} > {hmm_file} ')
	protHits = hmmsearchProtein(proteomeFile,hmm_file,sres_folder,hmm_evalue)
	if len(protHits) > 0:
		validated_clusters,validated_clusters_UIDS,protHits = find_clusters(protHits,nprot,prots_between,max_alien_prots,min_target_prots,min_cluster_coverage)
		os.remove(hmm_file)
		return(validated_clusters,validated_clusters_UIDS,sgenbankDict, description,proteomeFile,protHits)
	else:
		os.remove(hmm_file)
		return([],[],sgenbankDict, description,proteomeFile,protHits)

#blast only mode
def scan_cluster_Blastp(cluster_faa_file, nprot, prots_between,max_alien_prots, sgenbankFile, res_folder,evalue,qcov,scov, min_target_prots,min_cluster_coverage):
	#extract proteome
	if not os.path.exists(sgenbankFile):
		write_to_log(log_file,"Error: The genome genbank file could not be found.")
		sys.exit("Error: The genome genbank file could not be found.")
	try:
		sgenome_name = os.path.splitext(os.path.basename(sgenbankFile))[0]
		sres_folder = os.path.join(res_folder,sgenome_name)
		if not os.path.exists(sres_folder):
			os.mkdir(sres_folder)
		sgenbankDict, description = genbankToDict(sgenbankFile)
		proteomeFile = extractProteome(sgenbankDict,sres_folder,sgenome_name)
	except:
		write_to_log(log_file,f"Error: Is {sgenome_name} genome in genbank format?")
		sys.exit(f"Error: Is {sgenome_name} genome in genbank format?")
	#local Blast cluster_faa_file vs proteomeFile.
	blast_out = os.path.join(sres_folder,'blast.res')
	blast_command = f"blastp -subject {proteomeFile} -query {cluster_faa_file} -out {blast_out} -max_target_seqs {max_target} -word_size 3 -evalue {evalue} -outfmt '6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore sseq' 1>> {log_file}"
	write_to_log(log_file,f'Running blastp to retrieve homologs for cluster identification ({sgenbankFile}).')
	print(f'> Running blastp to retrieve homologs for cluster identification ({sgenbankFile}).')
	os.system(blast_command)
	if os.path.getsize(blast_out) == 0:
		write_to_log(log_file,'Error: Blast error, 0 bits result file...')
		print('Error: Blast error, 0 bits result file...')
		return('','','','','','')
	else:
		#read blast resutls
		blast_fields = ['query name', 'target name', 'pident', 'length', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'evalue', 'bitscore', 'sseq']
		blast_tabl = pd.read_table(blast_out, names = blast_fields)
		blast_tabl['qcov'] = [ (row['qend']-row['qstart'])/row['qlen'] for _,row in blast_tabl.iterrows()]
		blast_tabl['scov'] = [ (row['send']-row['sstart'])/row['slen'] for _,row in blast_tabl.iterrows()]
		blast_tabl = blast_tabl[(blast_tabl['qcov'] >= qcov ) & (blast_tabl['scov'] >= scov )]
		#read blast results for each protein
		blast_tabl.sort_values(['evalue'])
		#validate clusters
		validated_clusters,validated_clusters_UIDS,blast_tabl = find_clusters(blast_tabl,nprot,prots_between,max_alien_prots,min_target_prots,min_cluster_coverage)
		return(validated_clusters,validated_clusters_UIDS,sgenbankDict, description,proteomeFile,blast_tabl)

def extract_gb_region_for_clinker(sgenbankFile, sgenbankDict, validated_clusters_UIDS, gb_regions_for_clinker):
	sgenbankFile_basename = os.path.splitext(os.path.basename(sgenbankFile))[0]
	r_file_paths = []
	for cluster in validated_clusters_UIDS:
		start_coord = ''
		end_coord = ''
		for k in sgenbankDict:
			if cluster[0] in sgenbankDict[k].keys():
				cluster_replicon = k
				start_coord = min([sgenbankDict[k][uid][3] for uid in cluster])-50
				if start_coord < 0:
					start_coord = 0
				end_coord = max([sgenbankDict[k][uid][4] for uid in cluster])+50
				break
		if (start_coord == '') or (end_coord == ''):
			message = f'Error: Start or end extraction locations not found...'
			write_to_log(log_file, message)
			sys.exit(message)
		for record in SeqIO.parse(sgenbankFile,'gb'):
			if record.id == cluster_replicon:
				region_for_clinker = record[start_coord:end_coord]
				rname = f'{sgenbankFile_basename}__{cluster_replicon}__{start_coord}-{end_coord}.gb'
				rname_path = os.path.join(gb_regions_for_clinker,rname)
				r_file_paths.append(rname_path)
				SeqIO.write(region_for_clinker,rname_path,'gb')
	return(r_file_paths)

def align_target_proteins(qprot, target_for_alignment, proteins,alignments_folder):
	proteins_file = os.path.join(alignments_folder,re.sub('\\|','__',qprot)+'.fasta')
	alignment_file = os.path.join(alignments_folder,re.sub('\\|','__',qprot)+'_aln.fasta')
	if (os.path.exists(alignment_file)) and (os.path.getsize(alignment_file) > 0):
		message = f'--> MSA for {qprot} orthologs found...'
		print(message)
		write_to_log(log_file,message)
		return(alignment_file)
	else:
		message = f'--> Creating MSA for {qprot} orthologs...'
		print(message)
		write_to_log(log_file,message)
		to_fasta = []
		for protein in SeqIO.parse(proteins,'fasta'):
			pid = '|'.join(protein.id.split('|')[0:2])
			if pid in target_for_alignment:
				to_fasta.append(protein)
		SeqIO.write(to_fasta,proteins_file,'fasta')
		os.system(f'mafft --quiet --localpair --maxiterate 1000 {proteins_file} > {alignment_file}')
	return(alignment_file)

def get_id_dict_from_msa(qprot_alignment):
	sequences = {}
	i=0
	for r in SeqIO.parse(qprot_alignment,'fasta'):
		sequences[i] = (r.id,r.seq)
		i+=1
	id_dict = {}
	for i in range(0,len(sequences)):
		for j in range(i+1,len(sequences)):
			seq1_id,seq1 =  sequences[i]
			seq2_id,seq2 =  sequences[j]
			id_dict[(seq1_id,seq2_id)] = get_percent_id(seq1,seq2)
	return(id_dict)

def get_percent_id(seq1,seq2):
	id = 0
	alen = len(seq1)
	if alen == len(seq2):
		for i in range(alen):
			if (seq1[i] == seq2[i]) and (seq1[i] != '-') and (seq2[i] != '-'):
				id += 1
	else:
		sys.exit(f'Error: Sequences in the alignment have different length...')
	return(id/alen*100)

def get_alien_proteins_id_scores(cluster_dict, proteins, alignments_folder,evalue,id_dict):
	#identify all alien proteins
	alien_proteins = [ protein[0] for cluster in cluster_dict.values() for protein in cluster if protein[1] == ''] 
	alien_protein_records = []
	for protein in SeqIO.parse(proteins,'fasta'):
		if protein.id in alien_proteins:
			alien_protein_records.append(protein)
	if len(alien_protein_records) <= 1:
		return(id_dict)
	alien_protein_sequences = os.path.join(alignments_folder,'alien_proteins.fasta')
	blast_output = 'blast_res.tmp'
	SeqIO.write(alien_protein_records, alien_protein_sequences ,'fasta')
	blast_command = f"blastp -subject {alien_protein_sequences} -query {alien_protein_sequences} -out {blast_output} -word_size 3 -evalue {evalue} -outfmt '6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore sseq' 1>> {log_file}"
	print('\n> Running blastp to assign homologs between proteins found in the subject clusters but that were not present in the query cluster...')
	write_to_log(log_file,f'Running blastp to assign homologs between proteins found in the subject clusters but that were not present in the query cluster...')
	try:
		os.system(blast_command)
	except:
		message = f'Error: Blast search failed. {blast_command}'
		write_to_log(log_file,message)
		print(message)
	#Check blast results
	if os.path.getsize(blast_output) == 0:
		write_to_log(log_file,'Error: Remote Blast produced an error... Please try again or run in blast mode without HMM generation...')
		print('Error: Remote Blast produced an error... Please try again or run in blast mode without HMM generation...')
	#read blast resutls
	blast_fields = ['qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'evalue', 'bitscore', 'sseq']
	blast_tabl = pd.read_table(blast_output, names = blast_fields)
	blast_tabl = blast_tabl[ blast_tabl['qseqid'] != blast_tabl['sseqid'] ]
	blast_tabl['qcov'] = [ (row['qend']-row['qstart'])/row['qlen'] for _,row in blast_tabl.iterrows()]
	blast_tabl['scov'] = [ (row['send']-row['sstart'])/row['slen'] for _,row in blast_tabl.iterrows()]
	blast_tabl = blast_tabl[(blast_tabl['qcov'] >= qcov ) & (blast_tabl['scov'] >= scov )]
	os.remove(blast_output)
	processed =set()
	aprots = set(blast_tabl['qseqid'])
	for aprot in aprots:
		if not aprot in processed:
			target_for_alignment = list(blast_tabl[blast_tabl['qseqid'] == aprot]['sseqid'])
			target_for_alignment.append(aprot)
			processed = processed.union(set(target_for_alignment))
			alien_aignment_file = align_target_proteins(aprot, target_for_alignment, proteins,alignments_folder)
			id_dict.update(get_id_dict_from_msa(alien_aignment_file))
	return(id_dict)

#DP alignment works, but to improve
def compareByDP(cluster_i,cluster_j,id_dict, gap, mismatch):
	if cluster_i != cluster_j:
		nA = len(cluster_i)
		nB = len(cluster_j)
		max_len = max(nA,nB)
		dictA = {}
		dictB = {}
		dpMat = {}
		#gap penalty is substracted
		gap = abs(gap)
		mismatch = abs(mismatch)
		#index begin in 1 so that in the dpMat 0 indexes correspond to gap openings
		for i in range(1,nA+1):
			dictA[i] = cluster_i[i-1]
		for j in range(1,nB+1):
			dictB[j] = cluster_j[j-1]
		#initialize matrix
		for i in range(0,nA+1):
			for j in range(0,nB+1):
				if i == 0:
					dpMat[i,j] = [-gap*j,""]
				elif j == 0:
					dpMat[i,j] = [-gap*i,""]
				else:
					dpMat[i,j] = [0,""]
		#Fill dpMat
		#wrong orientation? half score?
		for i in range(1,nA+1):
			for j in range(1,nB+1):
				AGap = dpMat[i-1,j][0]-gap
				BGap = dpMat[i,j-1][0]-gap
				#orthologs?
				if ((dictA[i][0],dictB[j][0]) in id_dict.keys()) and (dictA[i][3] == dictB[j][3]):
					AB = dpMat[i-1,j-1][0]+id_dict[dictA[i][0],dictB[j][0]]
				elif (dictA[i][0],dictB[j][0]) in id_dict.keys():
					#wrong orientation
					AB = dpMat[i-1,j-1][0]+id_dict[dictA[i][0],dictB[j][0]]/2
				else:
					#not orthologs, prefere gap inserition instead of aligning different genes --> -20 
					AB = dpMat[i-1,j-1][0] - mismatch
				if (AB >= AGap) and (AB >= BGap):
					direction = "Match"
					dpMat[i,j] = [AB,direction]
				elif (AGap > AB) and (AGap >= BGap):
					direction = "AGap"
					dpMat[i,j] = [AGap,direction]
				elif (BGap > AB) and (BGap > AGap):
					direction = "BGap"
					dpMat[i,j] = [BGap,direction]
		#Trace back and generate alignment. Hacer una funcion que haga el traceback. Tiene que ir a la última posición dpMat[nA,nB] --> mirar la dirección, alinear seqs y ir a la siguiente direcciçon
		alignmentA = []
		alignmentB = []
		i = nA
		j = nB
		score = dpMat[i,j][0]
		while (i != 0) or (j != 0):
			if dpMat[i,j][1] == 'Match':
				alignmentA.append(dictA[i])
				alignmentB.append(dictB[j])
				i = i - 1
				j = j - 1
			elif dpMat[i,j][1] == 'AGap':
				alignmentA.append(dictA[i])
				alignmentB.append('----')
				i = i - 1
			elif dpMat[i,j][1] == 'BGap':
				alignmentA.append('----')
				alignmentB.append(dictB[j])
				j = j - 1
			elif (i == 0) and (dpMat[i,j][1] == ''):
				alignmentA.append('----')
				alignmentB.append(dictB[j])
				j = j - 1
			elif (j == 0) and (dpMat[i,j][1] == ''):
				alignmentA.append(dictA[i])
				alignmentB.append('----')
				i = i - 1
		alignmentA.reverse()
		alignmentB.reverse()
		return(score/max_len,alignmentA,alignmentB)
	else:
		return(100, cluster_i, cluster_j)

def reverse_cluster(cluster_j):
	cluster_j = [(a,b,c,-d) for a,b,c,d in cluster_j ]
	cluster_j.reverse()
	return(cluster_j)

def orient_clusters(ref_cluster,cluster_dict,id_dict,gap,mismatch):
	new_cluster_dict = {}
	for k in cluster_dict:
		cl = cluster_dict[k]
		if cl == ref_cluster:
			new_cluster_dict[k] = cl
		else:
			score, aln_cluster_i, aln_cluster_j = compareByDP(cl,ref_cluster,id_dict,gap,mismatch)
			rcl = reverse_cluster(cl)
			rscore, raln_cluster_i, raln_cluster_j = compareByDP(rcl,ref_cluster,id_dict,gap,mismatch)
			if rscore > score:
				new_cluster_dict[k] = rcl
			else:
				new_cluster_dict[k] = cl
	return(new_cluster_dict)

#phylo functions
def put_clusters_in_tree_leaves(tree,reoriented_cluster_dict):
	for leave in tree.get_terminals():
		leave.names = [leave.name]
		leave.clusters = [v for k,v in reoriented_cluster_dict.items() if leave.name == os.path.basename(os.path.splitext(k)[0])]
		if len(leave.clusters) != 1:
			sys.exit('Error: some leaves have the same name...')
	return(tree)

def get_parent(tree, node):
	if tree != node:
		path_to_node = tree.get_path(node)
		if len(path_to_node) == 1:
			return(tree)
		else:
			return(path_to_node[-2])

#Progresive alignment
def compareByDP_cluster_profiles(list_cluster_i,list_cluster_j,id_dict,gap,mismatch):
	nseqsA = len(list_cluster_i)
	nseqsB = len(list_cluster_j)
	nA = len(list_cluster_i[0])
	nB = len(list_cluster_j[0])
	max_len = max(nA,nB)
	dictA = {}
	dictB = {}
	dpMat = {}
	#gap penalty is substracted
	gap = abs(gap)
	mismatch = abs(mismatch)
	#index begin in 1 so that in the dpMat 0 indexes correspond to gap openings
	for i in range(1,nA+1):
		dictA[i] = {}
		for ci, cluster_i in enumerate(list_cluster_i):
			dictA[i].update({ci:cluster_i[i-1]})
	for j in range(1,nB+1):
		dictB[j] = {}
		for cj, cluster_j in enumerate(list_cluster_j):
			dictB[j].update({cj:cluster_j[j-1]})
	#initialize matrix
	for i in range(0,nA+1):
		for j in range(0,nB+1):
			if i == 0:
				dpMat[i,j] = [-gap*j,""]
			elif j == 0:
				dpMat[i,j] = [-gap*i,""]
			else:
				dpMat[i,j] = [0,""]
	#Fill dpMat
	#wrong orientation: half score
	for i in range(1,nA+1):
		for j in range(1,nB+1):
			AGap = dpMat[i-1,j][0]-gap
			BGap = dpMat[i,j-1][0]-gap
			#orthologs? Sum of scores, nseqsA, nseqsB
			AB = dpMat[i-1,j-1][0]
			sum_of_scores = 0
			comparisons = 0
			for sA in range(nseqsA):
				for sB in range(nseqsB):
					comparisons += 1
					if ((dictA[i][sA][0],dictB[j][sB][0]) in id_dict.keys()) and (dictA[i][sA][3] == dictB[j][sB][3]):
						sum_of_scores += id_dict[dictA[i][sA][0],dictB[j][sB][0]]
					elif (dictA[i][sA][0],dictB[j][sB][0]) in id_dict.keys():
						#wrong orientation
						sum_of_scores += id_dict[dictA[i][sA][0],dictB[j][sB][0]]/2
					elif (dictA[i][sA][0] == '----' ) and (dictB[j][sB] == '----'):
						sum_of_scores += -gap
					else:
						#not orthologs, prefere gap inserition instead of aligning different genes --> -20 
						sum_of_scores += -mismatch
			AB += sum_of_scores/comparisons
			if (AB >= AGap) and (AB >= BGap):
				direction = "Match"
				dpMat[i,j] = [AB,direction]
			elif (AGap > AB) and (AGap >= BGap):
				direction = "AGap"
				dpMat[i,j] = [AGap,direction]
			elif (BGap > AB) and (BGap > AGap):
				direction = "BGap"
				dpMat[i,j] = [BGap,direction]
	#Trace back and generate alignment.From last postion in dpMat[nA,nB] --> look direction, align seqs and continue
	alignmentA = []
	alignmentB = []
	i = nA
	j = nB
	score = dpMat[i,j][0]
	while (i != 0) or (j != 0):
		if dpMat[i,j][1] == 'Match':
			alignmentA.append(dictA[i])
			alignmentB.append(dictB[j])
			i = i - 1
			j = j - 1
		elif dpMat[i,j][1] == 'AGap':
			alignmentA.append(dictA[i])
			alignmentB.append('----')
			i = i - 1
		elif dpMat[i,j][1] == 'BGap':
			alignmentA.append('----')
			alignmentB.append(dictB[j])
			j = j - 1
		elif (i == 0) and (dpMat[i,j][1] == ''):
			alignmentA.append('----')
			alignmentB.append(dictB[j])
			j = j - 1
		elif (j == 0) and (dpMat[i,j][1] == ''):
			alignmentA.append(dictA[i])
			alignmentB.append('----')
			i = i - 1
	alignmentA.reverse()
	alignmentB.reverse()
	#reconstruct clusters
	reconstructed_alignment_A = []
	reconstructed_alignment_B = []
	for sA in range(nseqsA):
		aliA_sA = [gen[sA] if gen != '----' else '----' for gen in alignmentA ]
		reconstructed_alignment_A.append(aliA_sA)
	for sB in range(nseqsB):
		aliA_sB = [gen[sB] if gen != '----' else '----' for gen in alignmentB ]
		reconstructed_alignment_B.append(aliA_sB)
	return(score/max_len,reconstructed_alignment_A,reconstructed_alignment_B)

def cluster_progressive_alignment(tree, id_dict,gap,mismatch):
	#get the closests leaves
	tree2 = deepcopy(tree)
	while len(tree2.get_terminals()) > 1:
		for t in tree2.get_terminals():
			parent = get_parent(tree2,t)
			if parent.is_preterminal():
				siblings = parent.get_terminals()
				align_profile_sum_of_scores(tree2,siblings,id_dict,gap,mismatch)
				break
	return(tree2)

def align_profile_sum_of_scores(tree2,siblings,id_dict,gap,mismatch):
	if len(siblings) != 2:
		sys.exit('Error: unresolved tree')
	else:
		list_cluster_i = siblings[0].clusters
		list_cluster_j = siblings[1].clusters
		list_names_i = siblings[0].names
		list_names_j = siblings[1].names
	score, aln_cluster_i_list, aln_cluster_j_list = compareByDP_cluster_profiles(list_cluster_i,list_cluster_j,id_dict,gap,mismatch)
	parent = get_parent(tree2, siblings[0])
	parent.clusters = aln_cluster_i_list + aln_cluster_j_list
	parent.names = list_names_i + list_names_j
	for sibling in siblings:
		parent.collapse(sibling.name)
	pass

def put_alignment_in_tree(tree,m_cluster_alignment,m_cluster_alignment_names):
	for l in tree.get_terminals():
		name = l.name
		aln = [ aln_clust for aln_clust,aln_name in zip(m_cluster_alignment,m_cluster_alignment_names) if name == aln_name][0]
		l.alignment = aln 
	pass

#ITOL annotation
def get_itol_len(tree):
	a0 = tree.get_terminals()[0]
	aln_len = len(a0.alignment)
	max_sizes = [0]*aln_len
	spacer = 101*(aln_len-1)
	for l in tree.get_terminals():
		for gp in range(aln_len):
			if l.alignment[gp] != '----':
				gene_size = l.alignment[gp][2]
				if gene_size > max_sizes[gp]:
					max_sizes[gp] = gene_size
	region_len = sum(max_sizes) + spacer
	return(region_len,max_sizes,aln_len)

def generate_file_for_itol(tree,res_folder):
	shapes = {-1:'TL', 0:'RE', 1:'TR'}
	region_len,max_sizes,aln_len = get_itol_len(tree)
	print(region_len)
	itol_file = os.path.join(res_folder,'iTOLData.txt')
	with open(itol_file,'w+') as f:
		f.write("DATASET_DOMAINS\nSEPARATOR COMMA\nDATASET_LABEL, Multiple Cluster Alingment\nCOLOR,#ff0000\nDATA\n")
		#nodeID,protLen,SHAPE|START|END|COLOR|LABEL,
		for terminal in tree.get_terminals():
			name = terminal.name
			alignment = terminal.alignment
			start_coord = 1
			end = 0
			f.write(f'{name},{region_len}')
			for i in range(aln_len):
				if alignment[i] == '----':
					start_coord += max_sizes[i] + 101
				else:
					correction = max_sizes[i] - alignment[i][2]
					if alignment[i][3] == -1:
						end = start_coord + 49
						gene_shape = shapes[alignment[i][3]]
						label = alignment[i][0].split('|')[1]
						f.write(f',{gene_shape}|{start_coord}|{end}|#cb2c31|')
						start_coord = end + 1
						end = start_coord + alignment[i][2] - 50
						gene_shape = shapes[0]
						f.write(f',{gene_shape}|{start_coord}|{end}|#cb2c31|{label}')
						start_coord = end + 99	+ correction
					else:
						end = start_coord + alignment[i][2] - 50
						gene_shape = shapes[0]
						label = alignment[i][0].split('|')[1]
						f.write(f',{gene_shape}|{start_coord}|{end}|#008080|{label}')
						start_coord = end + 1
						end = start_coord + 49
						gene_shape = shapes[alignment[i][3]]
						f.write(f',{gene_shape}|{start_coord}|{end}|#008080|')	
						start_coord = end + 99 + correction
			f.write(f'\n')
	pass

def printClusters(reoriented_cluster_dict,clusterFile):
	with open(clusterFile,'w') as cf:
		for k,v in reoriented_cluster_dict.items():
			clusterText = " ".join( [ list(g)[0] for g in v]  )
			cf.write(f'{k}\t{clusterText}\n')


#######################################################################
#######################################################################
#accessory functions

def extract_cluster_dna(genbankFile,replicon_id,cstart,cend,res_folder):
	query_cluster_fna = os.path.join(res_folder,f'{replicon_id}_query_cluster.fna')
	for r in SeqIO.parse(genbankFile,'gb'):
		if (r.id == replicon_id) or (replicon_id == False):
			query_cluster = r[cstart:cend]
			SeqIO.write(query_cluster,query_cluster_fna,'fasta')
	return(query_cluster_fna)

def detect_error_in_tree(tree):
	tree2 = deepcopy(tree)
	while len(tree2.get_terminals()) > 1:
		closest_siblings = get_closest_leaves(tree2)
		if len(closest_siblings) != 2:
			return(f'ERROR',tree2)
		parent = get_parent(tree2, closest_siblings[0])
		for sibling in closest_siblings:
			parent.collapse(sibling.name)
	return(tree2)

def get_closest_leaves(tree):
	terminals = []
	for terminal in tree.get_terminals():
		terminals.append(terminal)
	distance = ''
	for i in range(len(terminals)):
		terminal = terminals[i]
		path = tree.get_path(terminal)
		for j in range(i+1, len(terminals)):
			target = list(tree.find_elements(terminals[j]))[0]
			dist = tree.distance(terminal,target)
			if distance == '':
				distance = dist
				if len(path) > 1:
					parent = tree.get_path(terminal)[-2]
				else:
					parent = tree
			elif dist < distance:
				distance = dist
				if len(path) > 1:
					parent = tree.get_path(terminal)[-2]
				else:
					parent = tree
	closest_siblings = parent.get_terminals()
	return(closest_siblings)

#######################################################################
#######################################################################

#Main
if __name__ == "__main__":
	args = parseArgs()
	printSoftName()
	skip_generate_HMM = False
	#result folder
	now = datetime.now()
	formated_time = now.strftime("%H%M%S_%d%m%Y")
	if args.res_folder:
		res_folder = args.res_folder
	else:
		res_folder = 'Results'
	print(f'The results will be saved in {res_folder} folder...')
	if not os.path.exists(res_folder):
		os.mkdir(res_folder)
	elif os.path.exists(res_folder) and args.overwrite:
		print(f'{res_folder} exists. Overwritting results...')
	else:
		res_folder = f'{res_folder}-{formated_time}'
		os.mkdir(res_folder)
	#
	print()
	#log file
	log_file = os.path.join(res_folder,'log.txt')
	log = open(log_file, 'w')
	log.write(f'{NAME} {VERSION}\n')
	log.flush()
	log.close()
	#File extension
	gbext = args.gbext
	#query
	if (args.qfile) and (not args.qcluster) and (not args.hmm_folder):
		genbankFile = args.qfile
		if not os.path.exists(genbankFile):
			write_to_log(log_file,'Error: Query genbank file not found...')
			sys.exit('Error: Query genbank file not found...')
		if args.repid:
			replicon_id = args.repid
		else:
			number_of_replicons = len(list(SeqIO.parse(genbankFile,'gb')))
			if not number_of_replicons == 1:
				write_to_log(log_file,f'Error: The ID of the replicon/contig containing the cluster is required... ')
				sys.exit(f'Error: The ID of the replicon/contig containing the cluster is required... ')
			else:
				replicon_id = False
		#cstart and cend
		if args.cstart:
			cstart = int(args.cstart)
		else:
			write_to_log(log_file,f'Error: A cluster start location is required...')
			sys.exit(f'Error: A cluster start location is required...')
		if args.cend:
			cend = int(args.cend)
		else:
			write_to_log(log_file,f'Error: A cluster end location is required...')
			sys.exit(f'Error: A cluster end location is required...')
		if cstart > cend:
			write_to_log(log_file,f'Error: Cluster end location smaller than cluster start location...')
			sys.exit(f'Error: Cluster end location smaller than cluster start location...')
		hmm_folder = os.path.join(res_folder,'HMM')
		if not os.path.exists(hmm_folder):
			os.mkdir(hmm_folder)
	elif (args.qcluster) and not (args.hmm_folder):
		genbankFile = args.qcluster
		if not os.path.exists(genbankFile):
			write_to_log(log_file,'Error: Query genbank file not found...')
			sys.exit('Error: Query genbank file not found...')
		replicon_id = False
		cstart = 0
		cend = ''
		hmm_folder = os.path.join(res_folder,'HMM')
		if not os.path.exists(hmm_folder):
			os.mkdir(hmm_folder)
	elif args.hmm_folder:
		genbankFile = False
		qcluster = False
		hmm_folder = args.hmm_folder
		skip_generate_HMM = True
	else:
		write_to_log(log_file,f'Error: A query cluster file in Genbank format is required...')
		sys.exit(f'Error: A query cluster file in Genbank format is required...')
	#
	#Subject genome
	if (args.sfile) and (not args.sfolder):
		if os.path.exists(args.sfile):
			sgenbankFiles = [args.sfile]
		else:
			write_to_log(log_file,'Error: Subject genbank file or folder containing subject genomes not found...')
			sys.exit('Error: Subject genbank file or folder containing subject genomes not found...')
	elif args.sfolder:
		if os.path.exists(args.sfolder):
			sfolder = args.sfolder
			sgenbankFiles = glob.glob(os.path.join(sfolder,f'*.{gbext}'))
			subject_genomes = len(sgenbankFiles)
			if subject_genomes == 0:
				write_to_log(log_file,'Error: Subject folder empty...')
				sys.exit('Error: Subject folder empty...')
		else:
			write_to_log(log_file,'Error: Subject folder not found...')
			sys.exit('Error: Subject folder not found...')		
	else:
		write_to_log(log_file,'Error: Subject genbank file or folder containing subject genomes not found...')
		sys.exit('Error: Subject genbank file or folder containing subject genomes not found...')
	
	
	#cluster alignment args...
	gap = int(args.gap)
	mismatch = int(args.mismatch)
	################################################################
	################################################################

	#blast arguments
	evalue = float(args.evalue)
	blastp_database = args.blastp_database
	#validate blast database
	if not blastp_database in BLAST_DBS:
		write_to_log(log_file,f'Error: {blastp_database} blastp database not available...')
		sys.exit(f'Error: {blastp_database} blastp database not available...')
	max_target = int(args.max_target)
	qcov = int(args.qcov)/100
	scov = int(args.scov)/100
	#hmmsearch arguments
	hmm_evalue = float(args.hmm_evalue)
	hmm_cover = int(args.hmm_cover)/100
	ali_cover = int(args.ali_cover)/100
	#cluster finder arguments
	#minimum of target genes found
	min_target_prots = int(args.min_target_prots)
	min_cluster_coverage = float(args.min_cluster_coverage)
	#
	################################################################
	################################################################
	print('\nStarting the analysis...')
	#Database for hmm generation: local
	local_blast_db = args.local_blast_db
	if not os.path.exists(f'{local_blast_db}.psq') and (local_blast_db != ''):
		write_to_log(log_file,f'Error: Local blast database ({local_blast_db}) not found')
		sys.exit(f'Error: Local blast database ({local_blast_db}) not found')
	elif os.path.isdir(local_blast_db):
		write_to_log(log_file,f'Error: Local blast database ({local_blast_db}) is a folder... It should be a DB name without the file extension (DB Prefix)...')
		sys.exit(f'Error: Local blast database ({local_blast_db}) is a folder... It should be a DB name without the file extension (DB Prefix)...')
	#if local database generated from subject genomes
	if (args.local_blast_db_subject) and (args.sfolder):
		#make blast db
		write_to_log(log_file,'Generating Blast database from genbank files...')
		print('> Generating Blast database from genbank files...')
		local_blast_db_path = os.path.join(res_folder,'blastDB')
		if not os.path.exists(local_blast_db_path):
			os.mkdir(local_blast_db_path)
		local_blast_db = os.path.join(local_blast_db_path,'blastDB')
		proteomes = []
		for subject in sgenbankFiles:
			sgenome_name = os.path.splitext(os.path.basename(subject))[0]
			genbankDict, description = genbankToDict(subject)
			proteomes.append(extractProteome(genbankDict,local_blast_db_path,sgenome_name))
		blastp_faa = os.path.join(local_blast_db_path,'blastp_db.faa')
		os.system(f'cat {" ".join(proteomes)} > {blastp_faa}')
		os.system(f'makeblastdb -dbtype prot -in {blastp_faa} -out {local_blast_db} 1>> {log_file}')	
	#
	print(f'-> Using local blast database: {local_blast_db}...')
	################################################################
	################################################################
	only_blastp = args.only_blastp
	#
	print(f'\nProcessing genomes...')
	subject_file_basenames = [ os.path.basename(f) for f in sgenbankFiles]
	if (genbankFile) and (not args.qcluster):
		ref_file_basename = os.path.basename(genbankFile) 
		if not ref_file_basename in subject_file_basenames:
			subject_file_basenames.append(genbankFile)
	elif args.qcluster and args.refgenome:
		ref_file_basename = os.path.basename(args.refgenome)
		if os.path.exists(args.refgenome):
			print(f'-> Using {args.refgenome} as reference genome...')
		else:
			print(f'-> {args.refgenome} not found- Using {sgenbankFiles[0]} as reference genome...')
			ref_file_basename = os.path.basename(sgenbankFiles[0])
	else:
		ref_file_basename = os.path.basename(sgenbankFiles[0])
	hit_table = pd.DataFrame(columns=['query name','target name', 'id_score','qcover','cluster','cluster_file'])
	gb_regions_for_clinker = os.path.join(res_folder,'GB_for_Clinker')
	if not os.path.exists(gb_regions_for_clinker):
		os.mkdir(gb_regions_for_clinker)
	if not only_blastp:
		if skip_generate_HMM:
			#hmm_folder
			HMMs = glob.glob(os.path.join(hmm_folder,'*.hmm'))
			nprot = len(HMMs)
			if nprot == 0:
				write_to_log(log_file,f'Error: No HMM file was found on {hmm_folder}...')
				sys.exit(f'Error: No HMM file was found on {hmm_folder}...')
		else:
			cluster_faa_file = extract_cluster(genbankFile,replicon_id,cstart,cend,res_folder)
			#generates HMMs for each protein in cluster. If local_blast_db='' it uses remote search (slow)
			HMMs = generate_HMM(cluster_faa_file,evalue,blastp_database,qcov,scov,res_folder,local_blast_db,hmm_folder)
			nprot = len(HMMs)
		#define the numbert of protein between consecutive genes
		if args.prots_between:
			prots_between = int(args.prots_between)
		else:
			prots_between = int(nprot/2)
		#max_alien_prots
		if args.max_alien_prots:
			max_alien_prots = int(args.max_alien_prots)
		else:
			max_alien_prots = int(nprot*3)
		#search for clusters
		if len(sgenbankFiles) == 1:
			print(f'Please use Only Blast method to run with a single genome...')
		for sgenbankFile in sgenbankFiles:
			print(f'\n-> Processing {sgenbankFile} ...')
			if AnnotatedGenome(sgenbankFile):
				validated_clusters,validated_clusters_UIDS,sgenbankDict, description,proteomeFile,protHits = scan_cluster_HMMER(HMMs, nprot, prots_between,max_alien_prots, sgenbankFile, res_folder, min_target_prots,min_cluster_coverage,hmm_evalue)
				if len(validated_clusters) == 0:
					print(f'--> No clusters found on {sgenbankFile}...')
					continue
				#get coordinates of genomic region to extract from genbank file --> for clinker
				region_file_path = extract_gb_region_for_clinker(sgenbankFile, sgenbankDict, validated_clusters_UIDS, gb_regions_for_clinker)
				#
				protHits.loc[:,'cluster_file'] = [ region_file_path[row['cluster']] for _,row in protHits.iterrows()]
				#protHits, use HMM acc as measure of prediction quality
				protHits = protHits.rename(columns={'acc': "id_score", "hmm cover": "qcover"}, errors="raise")
				columns_to_concat = ['query name','target name', 'id_score','qcover','cluster','cluster_file']
				protHits_subset = protHits.loc[:, columns_to_concat]
				# Identify rows where all values in the selected columns are NaN or empty
				rows_to_exclude = protHits_subset.isna().all(axis=1) | (protHits_subset == '').all(axis=1)
				# Filter protHits to exclude these rows
				protHits_filtered = protHits_subset[~rows_to_exclude]
				# Now concatenate with the filtered DataFrame
				if (not len(protHits_filtered) == 0) and (not len(hit_table) == 0):
					hit_table = pd.concat([hit_table, protHits_filtered])
				elif (len(hit_table) == 0) and (not len(protHits_filtered) == 0):
					hit_table = protHits_filtered
			else:
				print(f'--> {sgenbankFile} is not annotated....')
				continue
	else:
		cluster_faa_file = extract_cluster(genbankFile,replicon_id,cstart,cend,res_folder)
		nprot = len([1 for i in SeqIO.parse(cluster_faa_file,'fasta')])
		#define the numbert of protein between consecutive genes
		if args.prots_between:
			prots_between = int(args.prots_between)
		else:
			prots_between = int(nprot/2)
		#max_alien_prots
		if args.max_alien_prots:
			max_alien_prots = int(args.max_alien_prots)
		else:
			max_alien_prots = int(nprot/2)
		#search for clusters
		for sgenbankFile in sgenbankFiles:
			print(f'\n-> Processing {sgenbankFile} ...')
			validated_clusters,validated_clusters_UIDS,sgenbankDict, description,proteomeFile,blast_tabl = scan_cluster_Blastp(cluster_faa_file, nprot, prots_between,max_alien_prots, sgenbankFile, res_folder,evalue,qcov,scov, min_target_prots,min_cluster_coverage)
			if len(validated_clusters) == 0:
				print(f'--> No clusters found on {sgenbankFile}...')
				continue
			#get coordinates of genomic region to extract from genbank file --> for clinker
			region_file_path = extract_gb_region_for_clinker(sgenbankFile, sgenbankDict, validated_clusters_UIDS, gb_regions_for_clinker)
			#blast_tabl use pident as measure of
			blast_tabl.loc[:,'cluster_file'] = [ region_file_path[row['cluster']] for _,row in blast_tabl.iterrows()] 
			blast_tabl = blast_tabl.rename(columns={'pident': "id_score", "qcov": "qcover"}, errors="raise")
			if (len(hit_table) == 0) and (len(blast_tabl) > 0):
				hit_table = blast_tabl.loc[:,['query name','target name', 'id_score','qcover','cluster','cluster_file']]
			elif (not len(hit_table) == 0) and (len(blast_tabl) > 0):
				hit_table = pd.concat([hit_table, blast_tabl.loc[:,['query name','target name', 'id_score','qcover','cluster','cluster_file']]])
	#
	################################################################
	################################################################
	#analyze clusters
	######
	print('\n\nAnalyzing clusters...')
	write_to_log(log_file, 'Analyzing clusters...')
	if len(hit_table) == 0:
		sys.exit('No homologs were found.')
	hit_table.loc[:,'locus tag'] = [ row['target name'].split('|')[2] for _,row in hit_table.iterrows() ]
	hit_table.loc[:,'replicon_locus_tag'] = [ f"{row['target name'].split('|')[0]}|{row['target name'].split('|')[2]}" for _,row in hit_table.iterrows() ]
	hit_table = hit_table.reset_index(drop=True)
	region_files = glob.glob(os.path.join(gb_regions_for_clinker,'*.gb'))
	#here regions for clinker are read and compared based on the query proteins/cluster. The order and orientation of the genes are taken in account
	#function that returns the structure of the operon with the gene orientations
	#Function that makes the rc structure of the operon. 
	#cluster: contains locus_tag, query, len, strand, acc/%id vs query, qc
	# for proteins that are in the query but not in the subject --> ''
	alignments_folder = os.path.join(res_folder,'alignments')
	if not os.path.exists(alignments_folder):
		os.mkdir(alignments_folder)
	proteins = os.path.join(res_folder,'all_cluster_proteins.fasta')
	if os.path.exists(proteins):
		os.remove(proteins)
	#generate cluster dict
	ref_region = [rf if os.path.splitext(ref_file_basename)[0] in rf else region_files[0] for rf in region_files ][0]
	print(f'> The selected reference cluster was : {ref_region}')
	#get reference cluster
	ref_cluster = get_oriented_cluster_with_score(ref_region, hit_table,proteins)
	#cluster dictionary of all the found clusters
	cluster_dict = { region_file:get_oriented_cluster_with_score(region_file, hit_table,proteins) for region_file in region_files}
	##MAKE MSA for all proteins in cluster
	#For all cluster proteins a MSA is made and pairwise %ID is calcualted as a distance meassure to use in DP
	query_proteins = set(hit_table['query name'])
	id_dict = {}
	for qprot in query_proteins:
		qprot_hit_table = hit_table[hit_table['query name'] == qprot]
		#align and get %ID/Sim/Score
		target_for_alignment = set(qprot_hit_table['replicon_locus_tag'])
		qprot_alignment = align_target_proteins(qprot, target_for_alignment, proteins,alignments_folder)
		#get %id matrix for the MSA
		id_dict.update(get_id_dict_from_msa(qprot_alignment))
	#id_dict contains %id for all the comparisons, but the alien cluster proteins 
	#are not taken in account yet. 
	id_dict = get_alien_proteins_id_scores(cluster_dict, proteins, alignments_folder,evalue,id_dict)
	#generate reciprocal ids in id_dict
	id_dict.update({(k[1],k[0]):v for k,v in id_dict.items()})
	#id_dict contains the %ID off all ortholog comparisons to use as score in DP
	#align and clustering by DP
	#put all clusters with reference orientation
	reoriented_cluster_dict = orient_clusters(ref_cluster,cluster_dict,id_dict,gap,mismatch)
	#reverse? cluster_dict has '' if gene not found in the query sequence
	cluster_dict_keys = list(reoriented_cluster_dict.keys())
	pretty_cluster_dict_keys = [os.path.basename(os.path.splitext(k)[0]) for k in cluster_dict_keys]
	#
	clusterFile = os.path.join(res_folder,'clusters.txt')
	printClusters(reoriented_cluster_dict,clusterFile)
	#
	if len(reoriented_cluster_dict) > 200:
		sys.setrecursionlimit(len(reoriented_cluster_dict)*100)
	if len(reoriented_cluster_dict) > 1:
		#Generating distance cluster matrix!
		print('\nGenerating cluster distance tree...')
		print('> Generating distance cluster matrix...')
		idsMAT = {}
		dmatrix = []
		aln_matrix = {}
		#max score for normalization == 100
		for i in range(len(cluster_dict_keys)):
			matrixline = []
			for j in range(0,i+1):
				cluster_i = reoriented_cluster_dict[cluster_dict_keys[i]]
				cluster_j = reoriented_cluster_dict[cluster_dict_keys[j]]
				score, aln_cluster_i, aln_cluster_j = compareByDP(cluster_i,cluster_j,id_dict,gap,mismatch)
				dist = 1 - score/100
				#pretty_id_tags
				idsMAT[i,j] = pretty_cluster_dict_keys[i],pretty_cluster_dict_keys[j]
				aln_matrix[i,j] = aln_cluster_i, aln_cluster_j
				matrixline.append(dist)
			dmatrix.append(matrixline)
		##
		#hacer...
		print("> Constructing upgma tree...")
		write_to_log(log_file, 'Constructing upgma tree...')
		dm = TreeConstruction.DistanceMatrix(names=pretty_cluster_dict_keys, matrix= dmatrix )
		constructor = TreeConstruction.DistanceTreeConstructor()
		tree = constructor.upgma(dm)
		#tree = constructor.nj(dm) #for 3 sequences gives unresolved tree
		treeFile = os.path.join(res_folder,"tree.nhx")
		print(f'-> The tree file is located in: {treeFile}')
		write_to_log(log_file, f'The tree file is located in: {treeFile}')
		Phylo.write(tree, treeFile, "nexus")
		#
		print(f'\n> Generating Multiple cluster alignment with global Dynamic programming')
		tree = put_clusters_in_tree_leaves(tree,reoriented_cluster_dict)
		tree2 = cluster_progressive_alignment(tree, id_dict,gap,mismatch)
		m_cluster_alignment = tree2.clusters
		m_cluster_alignment_names = tree2.names
		put_alignment_in_tree(tree,m_cluster_alignment,m_cluster_alignment_names)
		print(f'> Generating Multiple cluster alignment ITOL annotation for the cluster distance tree: {os.path.join(res_folder,"iTOLData.txt")}')
		#	
		generate_file_for_itol(tree,res_folder)
		write_to_log(log_file, f'Generating Multiple cluster alignment ITOL annotation for the cluster distance tree: {os.path.join(res_folder,"iTOLData.txt")}')
	else:
		print(f'Only one cluster was found... No tree will be generated...')
	print(f'\nAnalysis done! Please find your results on the {res_folder} folder...\nThank you for using {NAME}!\n\n')