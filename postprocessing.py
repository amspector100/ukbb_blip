"""
Combines the MCMC results with the SuSiE results
and creates the final outputs except for the plots.
"""
import copy
import json
import numpy as np
import pandas as pd

import sys
import os
import parser
import time

import preprocessing
from preprocessing import elapsed
import getdata

ERRORS = ['fdr']
ALL_TRAITS = [
	'body_HEIGHTz', 
	'disease_CARDIOVASCULAR', 
	'biochemistry_HDLcholesterol',
	'biochemistry_LDLdirect'
]
METHOD_NAMES = {
	"susie_blip":"blip + susie",
	"susie":"susie",
	"susie-indiv-only":"susie (no groups)",
}

def create_final_json(trait, time0, args):

	## Load all rejections in json format
	# Susie, susie + blip, susie indiv only
	with open(f"output/polyfun_wrapper/rejections/{trait}.json", 'r') as thefile:
		all_rej = json.load(thefile)
	
	## Task 1: get other metadata by merging properly
	# Step 1a: find set of all snps
	all_snps = set()    
	methods = sorted(list(all_rej.keys()))
	for method in methods:
		rej = all_rej[method]
		# Find list of all snps
		for x in all_rej[method]:
			all_snps = all_snps.union(set(x['snps']))
	all_snps = sorted(list(all_snps))

	# Step 1b: Form dataframe of metadata
	meta_df_fname = f"output/final/rejections/metadata_{trait}.csv"
	if os.path.exists(meta_df_fname) and not args.get("regen_metadata", [False])[0]:
		print("Not regenerating final metadata. This saves time but may cause errors.")
		print("If there are errors, try adding the argument '--regen_metadata True'.")
		final_meta = pd.read_csv(meta_df_fname, index_col=0)
	else:
		print("Regenerating metadata. This may require downloading large files and it may take a while.")
		meta_df = pd.DataFrame(index=all_snps)
		meta_df['snp'] = meta_df.index
		split = meta_df['snp'].str.split(".", expand=True)
		meta_df['CHR'] = split[0].astype(int)
		meta_df['BP'] = split[1].astype(int)
		meta_df['ALLELE0'] = split[2]
		meta_df['ALLELE1'] = split[3]
		meta_df = meta_df.drop('snp', axis='columns')
		# Step 1c: merge with other metadata, by chromosome
		all_meta = []
		
		# helpers to load sumstats quickly
		getdata.pull_main_data(download_bolt=True, download_susie=False)
		sumstats_file = f"data/polyfun_results/bolt_337K_unrelStringentBrit_MAF0.001_v3.{trait}.bgen.stats.gz"
		chromepos_file = f"main_cache/chromepos/{trait}_chromepos.csv"
		chrome_starts = pd.read_csv(chromepos_file)
		chrome_starts = chrome_starts.set_index("CHR")
		chrome_starts.index = chrome_starts.index.astype(int)
		
		# loop through chromosomes
		for chrome in range(1, 23):
			print(f"Merging with metadata at chrome={chrome}, time={elapsed(time0)}.")
			subset = meta_df.loc[meta_df['CHR'] == chrome].copy()
			orig_n = subset.shape[0]
			sumstats = pd.read_csv(
				sumstats_file,
				sep="\t",
				header=0,
				usecols=[0,1,2,3,4,5,6,7],
				nrows=chrome_starts.loc[chrome, 'counts'],
				skiprows=chrome_starts.loc[chrome, 'starts']
			)
			sumstats.columns = ['SNP', 'CHR', 'BP', 'GENPOS', 'ALLELE1', 'ALLELE0', 'A1FREQ', 'INFO']
			print(f"Finished loading sumstats for chr={chrome} at {elapsed(time0)}.")
			subset = pd.merge(
				subset.reset_index(), sumstats, on=['CHR', 'BP', 'ALLELE0', 'ALLELE1'], how='left', copy=False
			).set_index('index')
			print(f"Finished merging with sumstats for chr={chrome} at {elapsed(time0)}.")
			new_n = subset.shape[0]
			if new_n != orig_n:
				raise ValueError(f"After merging chrome={chrome}, trait={trait}, subset shape changed from {new_n} to {orig_n}.")        
			# Append
			all_meta.append(subset)

		# Concatenate, add metadata, and save
		final_meta = pd.concat(all_meta, axis='index')
		final_meta.to_csv(meta_df_fname)

	# Save final dictionary with metadata
	final_dict = {}
	final_meta_snps = set(list(final_meta.index))
	for method in methods:
		final_dict[method] = []
		for x in all_rej[method]:
			snps = x['snps']
			#print(method)
			#print(set(snps) - set(snps).intersection(final_meta_snps))
			final_dict[method].append(
				dict(
					pep=x['pep'],
					SNP_IDS=snps,
					CHR=final_meta.loc[snps, 'CHR'].tolist(),
					BP=final_meta.loc[snps, 'BP'].tolist(),
					SNPS=final_meta.loc[snps, 'SNP'].tolist(),
					ALLELE0=final_meta.loc[snps, 'ALLELE0'].tolist(),
					ALLELE1=final_meta.loc[snps, 'ALLELE1'].tolist(),
					A1FREQ=final_meta.loc[snps, 'A1FREQ'].tolist(),
				)
			)
	with open(f"output/final/rejections/rejections_{trait}.json", "w") as savefile:
		savefile.write(json.dumps(final_dict))

def main(args):

	### Main arguments:
	# 1. create_final_json: whether or not to add metadata to raw outputs
	# 2. traits (list of traits)
	time0 = time.time()
	args = parser.parse_args(args)
	traits = args.get("traits", ALL_TRAITS)
	# Add all metadata, combine all outputs
	if args.get("create_final_json", [True])[0]:
		for trait in traits:
			create_final_json(trait, time0, args)

	# Create final plots
	#main_plots(traits, time0, args)

if __name__ == '__main__':
	main(sys.argv)