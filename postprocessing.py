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

# Plotting
import matplotlib.pyplot as plt
from plotnine import *


ERRORS = ['fdr']
ALL_TRAITS = [
	'body_HEIGHTz', 
	'disease_CARDIOVASCULAR', 
	'biochemistry_HDLcholesterol',
	'biochemistry_LDLdirect'
]
COLORS = [
	# 'cornflowerblue',
	# 'mediumaquamarine', #'aquamarine is pretty but hard to see, similarly for springgren
	# 'orangered',
	# 'pink',
	# 'gray',
	# 'gold'
	'#019E73', # blue-ish green
	'cornflowerblue', # blue
	"orangered", #'#D55E01', # red/vermillon
	'#F0E442', # yellow
	'#BBBBBB', # gray
	'#000000', # black
]
METHOD_NAMES = {
	"susie_blip":"blip + susie",
	"susie":"susie",
	"susie-indiv-only":"susie (no groups)",
}

def create_final_json(trait, time0):

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
	if False:#os.path.exists(meta_df_fname):
		final_meta = pd.read_csv(meta_df_fname, index_col=0)
	else:
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

def main_plots(traits, time0, args):

	# Read final data
	loci_df_list = []
	for trait in traits:
		with open(f"output/final/rejections/rejections_{trait}.json", "r") as file:
			rejections = json.load(file)

		# Load metadata about positions in chromosomes
		chromepos_file = f"main_cache/chromepos/{trait}_chromepos.csv"
		chrome_starts = pd.read_csv(chromepos_file)
		chrome_starts = chrome_starts.set_index("CHR")
		chrome_starts.index = chrome_starts.index.astype(int)
		chrome_starts['cum_bp'] = np.cumsum(chrome_starts['max_bp'])
		chrome_starts['cum_bp'] = chrome_starts['cum_bp'] - chrome_starts['max_bp']

		# Turn into dataframe
		snp_list = []
		loci_list = []
		methods = list(rejections.keys())
		for method in methods:
			for x in rejections[method]:
				pep = x['pep']
				gsize = len(x['SNP_IDS'])
				chrome = x['CHR'][0]
				bp = x['BP'][0] + chrome_starts.loc[chrome, 'cum_bp']
				loci = int(bp / int(1e6) - 0.5)
				loci_list.append(
					[method, chrome, loci, 1, 1 / gsize, gsize, pep]
				)
				for j in range(len(x['BP'])):
					snp_list.append(
						[method, x['BP'][j], 1 / gsize, gsize, pep]
					)

		## Aggregated plots by loci for each trait
		loci_df = pd.DataFrame(
			loci_list, columns=['method', 'CHR', 'loci', 'ndisc', 'gap', 'gsize', 'pep']
		)
		loci_df['Method'] = loci_df['method'].map(METHOD_NAMES)
		base_title = f"Discoveries for {trait.split('_')[1]}"
		for loci_subset, point_size, point_stroke, col_size, title, save_postfix in zip(
			[loci_df, loci_df.loc[loci_df['CHR'] == 1]],
			[0.05, 0.8],
			[0.05, 0.25],
			[5, 0.5],
			[base_title, base_title+" for Chromosome 1"],
			['', '_chrome1']

		):
			chrome1 = save_postfix != ''
			g = (
				ggplot(
					loci_subset.loc[
						(loci_subset['method'] != 'susie-indiv-only') &
						(~loci_subset['method'].str.contains('lfdr'))
					], 
					aes(
						x='loci',
						y='gap',
						fill='Method'
					)
				) +
				facet_wrap("Method", nrow=4) +
				geom_col(position='stack', size=col_size) + #, color='black', size=0.2) +
				labs(
					x="Region Start (mBP)",
					y="Group-Adjusted Num. Disc. per Locus",
					title=title,
				) +
				#geom_text(aes(x='ld'), label="'", y=-0.1,  size=8, color='black') +
				#scale_color_manual(values=COLORS) +
				scale_fill_manual(values=COLORS) 
			)
			if chrome1:
				g += geom_point(position='stack', color='black', size=point_size, stroke=point_stroke)
			g.save(f"output/final/plots/{trait}_ld_plot{save_postfix}.png", dpi=500)

		# Append
		loci_df['trait'] = trait.split("_")[1]
		loci_df_list.append(loci_df)

		# ## Plot for each discovered group / snp
		# snp_df = pd.DataFrame(
		# 	snp_list, columns=['method', 'bp', 'gap', 'group_size', 'pep']
		# )
	
	# Concatenate and plot group sizes
	MAX_SIZE = 25
	breaks = [0, 1, 2, 3, 5, 8, 11, 15, 20, MAX_SIZE]
	all_loci_df = pd.concat(loci_df_list, axis='index')
	all_loci_df['bin'] = pd.cut(all_loci_df['gsize'], bins=breaks, right=True)
	if not args.get("include_lfdr", [0])[0]:
		all_loci_df = all_loci_df.loc[~all_loci_df['Method'].str.contains("lfdr")]

	gsize_df = pd.DataFrame(all_loci_df.groupby(['bin', 'Method', 'trait'])['bin'].count())
	gsize_df = gsize_df.rename({'bin':'count'}, axis='columns').reset_index()

	title = f"Discovered Group Sizes"
	g = (
		ggplot(gsize_df, aes(x='bin', y='count', fill='Method'))
			+ geom_col(position='dodge')
			+ labs(
				x="Group Size", y="Number of Rejections", 
				title=title
			)
		+ facet_wrap("~trait", scales="free_y")
		+ scale_fill_manual(values=COLORS)
		+ theme_bw()
	    + theme(axis_text_x=element_text(angle = 35), figure_size=(6,3))
		+ theme(axis_text_x=element_text(angle = 35))
	)
	g.save(f"output/final/plots/groupsizes.jpg", dpi=500)

	# Plot overall group-adjusted power
	gap_df = all_loci_df.groupby(['Method', 'trait'])['gap'].sum().reset_index()
	g = (
	    ggplot(gap_df, aes(x="Method", y="gap", fill="Method")) 
	    + geom_col(position='dodge') 
	    + facet_wrap("~trait", scales='free_y')
	    + scale_fill_manual(values=COLORS)
	    + theme_bw()
	    + theme(subplots_adjust={'wspace':0.15})
	    + theme(axis_text_x=element_text(angle = 35), figure_size=(6,3))
	    + labs(title=f"Resolution Adjusted Power on UK Biobank, N=337K", y="Res. Adj. Num. Disc.")
	)
	g.save(f"output/final/plots/gap.jpg", dpi=500)





def main(args):

	### Main arguments:
	# 1. create_final_json: whether or not to add metadata to raw outputs
	# 2. traits (list of traits)
	# 3. include_lfdr: include results from local fdr control

	time0 = time.time()
	args = parser.parse_args(args)
	traits = args.get("traits", ALL_TRAITS)
	# Add all metadata, combine all outputs
	if args.get("create_final_json", [True])[0]:
		for trait in traits:
			create_final_json(trait, time0)

	# Create final plots
	#main_plots(traits, time0, args)

if __name__ == '__main__':
	main(sys.argv)