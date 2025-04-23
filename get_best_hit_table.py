#!/usr/bin/env python3
#requirements Blast >2.13
NAME = "Scan cluster: Get Best Hits Table"
VERSION = "1.0"
REF = "\n   Not published"
GITHUB="https://github.com/maurijlozano"

#Imports
#Python modules
import argparse, sys, os, glob
from math import log10
import seaborn as sns
from matplotlib import pyplot as plt

#Third party modules
import pandas as pd

#Functions
def parseArgs():
	'''
	Argument parsing is done.
	'''
	parser = argparse.ArgumentParser(description='This uses the blast or hmmsearch results to produce a best hit table for each gene in the cluster.')
	#Genome with the query cluster, gb file
	parser.add_argument("-o", "--Results_folder",help="Results folder name.", dest="res_folder", action='store', default = 'Results')
	parser.add_argument("-g", "--Genomes_folder",help="Folder with the genomes used for Scan cluster analysis.", dest="g_folder", action='store', default = 'Genomes')
	parser.add_argument("--Genebank_Extension",help="Genebank file extension. Default 'gb'", dest="ext", action='store', default = 'gb')
	parser.add_argument("--HMM_Results",help="Set on if Scan cluster was ran using HMMsearch.", dest="hmm_mode", action='store_true')
	#Optional arguments for blast and HMM
	parser.add_argument("--hmm_evalue",help="E-value cut-off for hmmsearch. Default = 0.00001", dest="hmm_evalue", action='store', default=0.00001)
	parser.add_argument("--hmm_cover",help="HMM coverage cut-off for hmmsearch. Default = 45", dest="hmm_cover", action='store', default=45)
	parser.add_argument("--Blast_evalue",help="E-value cut-off for remote blastp, used to retrieve homologs for HMM generation.", dest="evalue", action='store', default=0.00001)
	parser.add_argument("--Blast_qcov",help="Query coverage percent for Blastp search. Default=45", dest="qcov", action='store', default=45)
	parser.add_argument("--Blast_scov",help="Subject coverage percent for Blastp search. Default=45", dest="scov", action='store', default=45)
	args = parser.parse_args()
	return args

def hmmsearchProtein(hmmf,hmm_evalue,hmm_cover):
	tableHeaders = ['target name', 'taccession', 'tlen', 'query name', 'qaccession', 'qlen', 'E-value', 'score', 'bias', '#', 'of', 'c-Evalue', 'i-Evalue', 'dom score', 'dom  bias', 'hmm from', 'hmm to', 'ali from', 'ali to', 'env from', 'env to', 'acc', 'description of target']
	hmmHits = pd.read_table(hmmf, names=tableHeaders, comment='#', sep='\\s+').sort_values(["target name","E-value"])
	hmmHits = hmmHits.sort_values("E-value")
	hmmHitsPerGene = hmmHits.groupby(['query name'], as_index=False, ).first().sort_values("E-value")
	hmmHitsPerGene = hmmHitsPerGene[hmmHitsPerGene['E-value'] < hmm_evalue]
	hmmHitsPerGene.loc[:,'hmm cover'] = [(row['hmm to']-row['hmm from'])/row['qlen'] for _,row in hmmHitsPerGene.iterrows()]
	hmmHitsPerGene = hmmHitsPerGene[hmmHitsPerGene['hmm cover'] >= hmm_cover]
	hmmHitsPerGene = hmmHitsPerGene.reset_index(drop=True)
	return(hmmHitsPerGene)

def hmmsearchProteinCount(hmmf,hmm_evalue,hmm_cover):
	tableHeaders = ['target name', 'taccession', 'tlen', 'query name', 'qaccession', 'qlen', 'E-value', 'score', 'bias', '#', 'of', 'c-Evalue', 'i-Evalue', 'dom score', 'dom  bias', 'hmm from', 'hmm to', 'ali from', 'ali to', 'env from', 'env to', 'acc', 'description of target']
	hmmHits = pd.read_table(hmmf, names=tableHeaders, comment='#', sep='\\s+').sort_values(["target name","E-value"])
	hmmHits = hmmHits.sort_values("E-value")
	hmmHits = hmmHits[hmmHits['E-value'] < hmm_evalue]
	hmmHits.loc[:,'hmm cover'] = [(row['hmm to']-row['hmm from'])/row['qlen'] for _,row in hmmHits.iterrows()]
	hmmHits = hmmHits[hmmHits['hmm cover'] >= hmm_cover]
	hmmHitsPerGene = hmmHits.groupby(['query name'], as_index=False, )['target name'].count()
	hmmHitsPerGene = hmmHitsPerGene.reset_index(drop=True)
	return(hmmHitsPerGene)


#Main
if __name__ == "__main__":
	args = parseArgs()
	res_folder = args.res_folder
	if not os.path.exists(res_folder):
		sys.exit(f'Results folder not found: {res_folder}')
	#
	g_folder = args.g_folder
	ext = args.ext
	if not os.path.exists(g_folder):
		sys.exit(f'Results folder not found: {g_folder}')
	genomes = glob.glob(f'{g_folder}/*.{ext}')
	genomesNames = [ os.path.splitext(os.path.basename(g))[0] for g in genomes] 
	#
	hmm_mode = args.hmm_mode
	if hmm_mode:
		print(f'Running for results of Scan cluster HMM mode...')
	else:
		print(f'Running for results of Scan cluster only blast mode...')
	#Subject genome
	################################################################
	################################################################
	#blast arguments
	evalue = float(args.evalue)
	qcov = int(args.qcov)/100
	scov = int(args.scov)/100
	#hmmsearch arguments
	hmm_evalue = float(args.hmm_evalue)
	hmm_cover = int(args.hmm_cover)/100
	################################################################
	################################################################
	if hmm_mode:
		tableHeaders = ['target name', 'taccession', 'tlen', 'query name', 'qaccession', 'qlen', 'E-value', 'score', 'bias', '#', 'of', 'c-Evalue', 'i-Evalue', 'dom score', 'dom  bias', 'hmm from', 'hmm to', 'ali from', 'ali to', 'env from', 'env to', 'acc', 'description of target']
		hmm_files = glob.glob(f'{res_folder}/*/hmmsearch.txt')
		#
		best_hits = pd.DataFrame(columns=tableHeaders)
		count_Hits = pd.DataFrame(columns=['query name','target name','Genome'])
		for hmmf in hmm_files:
			bH = hmmsearchProtein(hmmf,hmm_evalue,hmm_cover)
			bH.loc[:,'Genome'] = os.path.split(os.path.split(hmmf)[0])[1]
			countH = hmmsearchProteinCount(hmmf,hmm_evalue,hmm_cover)
			countH.loc[:,'Genome'] = os.path.split(os.path.split(hmmf)[0])[1]
			if len(best_hits) > 0 and len(bH) > 0:
				best_hits = pd.concat([best_hits,bH])
			elif len(bH) > 0:
				best_hits = bH
			if len(count_Hits) > 0 and len(countH) > 0:
				count_Hits = pd.concat([count_Hits,countH])
			elif len(countH) > 0:
				count_Hits = countH
		#
	else:
		blastRes = glob.glob(f'{res_folder}/*/blast.res')
		if len(blastRes) == 0:
			sys.exit(f'Error: blast results not found...')
		blast_fields = ['query name', 'target name', 'pident', 'length', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'evalue', 'score', 'sseq']
		best_hits = pd.DataFrame(columns=blast_fields)
		count_Hits = pd.DataFrame(columns=['query name','target name','Genome'])
		for blast_out in blastRes:
			blast_tabl = pd.read_table(blast_out, names = blast_fields)
			blast_tabl['qcov'] = [ (row['qend']-row['qstart'])/row['qlen'] for _,row in blast_tabl.iterrows()]
			blast_tabl['scov'] = [ (row['send']-row['sstart'])/row['slen'] for _,row in blast_tabl.iterrows()]
			blast_tabl = blast_tabl[(blast_tabl['qcov'] >= qcov ) & (blast_tabl['scov'] >= scov )]
			#read blast results for each protein
			blast_tabl.sort_values(['evalue'])
			blast_tabl_bh = blast_tabl.groupby(['query name'], as_index=False, ).first().sort_values("E-value")
			blast_tabl_count = blast_tabl.groupby(['query name'], as_index=False, )['target name'].count()
			blast_tabl_bh['Genome'] = os.path.split(os.path.split(blast_out)[0])[1]
			blast_tabl_count['Genome'] = os.path.split(os.path.split(blast_out)[0])[1]
			blast_tabl_bh = blast_tabl_bh.reset_index(drop=True)
			blast_tabl_count = blast_tabl_count.reset_index(drop=True)
			#
			if len(best_hits) > 0 and len(blast_tabl_bh) > 0:
				best_hits = pd.concat([best_hits,blast_tabl_bh])
			elif len(blast_tabl_bh) > 0:
				best_hits = blast_tabl_bh
			if len(count_Hits) > 0 and len(blast_tabl_count) > 0:
				count_Hits = pd.concat([count_Hits,blast_tabl_count])
			elif len(blast_tabl_count) > 0:
				count_Hits = blast_tabl_count
	####
	best_HITS = best_hits.groupby(['query name','Genome'], as_index=False).first()
	best_HITS.to_csv(f'{res_folder}_best_hits.csv')
	count_Hits.to_csv(f'{res_folder}_count_hits.csv')
	#presence absence matrix
	#best_HITS['LogE'] = [ log10(e) if e >= 1e-100 else 1e-100 for e in best_HITS['E-value'] ]
	wideTable = best_HITS.pivot(index='Genome', columns='query name', values='score')
	wideTable = wideTable.fillna(0)
	missingGenomes = [ g for g in wideTable.index if not g in genomesNames] 
	if len(missingGenomes) > 0:
		missingGenomesWideTable = pd.DataFrame(0, index=missingGenomes, columns=wideTable.columns)
		wideTable = pd.concat([wideTable,missingGenomesWideTable])
	wideTable.to_csv(f'{res_folder}_best_hits_wide.csv')
	fig_width = len(wideTable.columns)*2
	fig_heigh = len(wideTable.index)/2
	sns.set_theme(font_scale=2)
	g = sns.clustermap(wideTable, figsize=(fig_width,fig_heigh), method='average', metric='euclidean', standard_scale=1, row_cluster=True, col_cluster=True, cmap="Blues", dendrogram_ratio=0.07)
	plt.savefig(f'{res_folder}_best_hits_wide.svg')
	#plt.show()
	plt.close()

	wideTable_count = count_Hits.pivot(index='Genome', columns='query name', values='target name')
	wideTable_count = wideTable_count.fillna(0)
	missingGenomes = [ g for g in wideTable_count.index if not g in genomesNames] 
	if len(missingGenomes) > 0:
		missingGenomeswideTable_count = pd.DataFrame(0, index=missingGenomes, columns=wideTable_count.columns)
		wideTable_count = pd.concat([wideTable_count,missingGenomeswideTable_count])
	wideTable_count.to_csv(f'{res_folder}_count_Hits_wide.csv')
	fig_width = len(wideTable_count.columns)*2
	fig_heigh = len(wideTable_count.index)/2
	sns.set_theme(font_scale=2)
	g = sns.clustermap(wideTable_count, figsize=(fig_width,fig_heigh), method='average', metric='euclidean', standard_scale=None, row_cluster=True, col_cluster=True, cmap="Blues", dendrogram_ratio=0.07)
	plt.savefig(f'{res_folder}_count_Hits_wide.svg')
	#plt.show()
	plt.close()

print('Done!')