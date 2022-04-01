import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse
import time
from tqdm import tqdm
import pyblip
from pyblip.create_groups import CandidateGroup

import copy

from preprocessing import elapsed

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns

# Saving output
import json

NUM_COLS = [f'ALPHA{i}' for i in range(1, 11)]
COLORS = ['cornflowerblue', 'orangered', 'forestgreen']

def load_alphas(outcome, funct=False, time0=None):
	if funct is False:
		postfix = 'nonfunct'
	else:
		postfix = 'funct'

	filename = f'data/polyfun_results/alphas/{outcome}.{postfix}.alphas.parquet'
	alphas = pd.read_parquet(filename, engine='fastparquet')
	print(f"Finished loading alphas for {outcome} at {elapsed(time0)}.")
	alphas_np = alphas[NUM_COLS].values
	return alphas, alphas_np

def fast_sequential_peps(alphas, q=0.05, thresh=0.25, max_size=25):
	"""
	Computes sequential peps of Alphas, filtering out
	all that are over thresh.
	"""

	p = alphas.shape[0]
	L = alphas.shape[1]
	cumalphas = np.zeros((p + 1,L))
	cumalphas[1:(p+1)] = np.cumsum(alphas, axis=0)

	# Compute successive groups of size m
	all_PEPs = {}
	for m in tqdm(list(range(max_size))):
		log_cumdiffs = np.log(
			1 - (cumalphas[(m+1):(p+1)] - cumalphas[:int(p-m)])
		)
		PEPs = np.exp(log_cumdiffs.sum(axis=1))
		all_PEPs[m] = PEPs
		
	# the index is the first (smallest) variable in the group which has size m
	active_inds = {}
	for m in range(max_size):
		active_inds[m] = np.where(all_PEPs[m] < thresh)[0]

	# This iteratively updates the list elim_inds so that
	# when consider the set of groups of size m+1, 
	# elim_inds are all the indices that are redundant
	elim_inds = set(np.where(all_PEPs[0] < q)[0].tolist())
	for m in range(1,max_size):
		# If index j is eliminated for level m-1, indexes j and j-1
		# are eliminated for index m
		elim_inds = elim_inds.union(set([x-1 for x in elim_inds]))
		# Update active_inds[m]
		update = set(active_inds[m].tolist()) - elim_inds
		active_inds[m] = np.array(list(update))
		# At this level, account for groups with PEP < q
		elim_inds = elim_inds.union(set(
			np.where(all_PEPs[m] < q)[0].tolist()
		))
		#print(f"After filtering, for m={m}, size={active_inds[m].shape}")


		# diffs = active_inds[m].reshape(-1, 1) - single_inds.reshape(1, -1)
		# diffs = np.abs(diffs).min(axis=1)
		# inclusion = diffs > m
		# active_inds[m] = active_inds[m][inclusion]

	return all_PEPs, active_inds

def blip_results(
		alphas,
		alphas_index,
		single_peps,
		q, 
		rej_susie_meta, # susie discoveries
		thresh, 
		rel_thresh,
		max_size,
		time0
	):

	# Step 1: Only include relevant snps with marginal pep < rel_thresh
	p = alphas.shape[0]
	L = alphas.shape[1]
	rel_inds = single_peps <= rel_thresh
	nrel = np.sum(rel_inds)
	relevant = alphas[rel_inds]
	arangep = np.arange(p)

	# Helper function which maps original index to new index
	def orig2new(ind):
		if rel_inds[ind]:
			return int(np.sum(rel_inds[0:ind]))
		else:
			return int(nrel + ind)

	# Step 2: compute sequential PIPs
	cand_groups = pyblip.create_groups.sequential_groups(
		susie_alphas=relevant.T, q=q, max_pep=thresh, max_size=max_size
	)
	print(f"Finished computing sequential cand_groups at time={elapsed(time0)}.")

	# Step 3: Add SuSiE's groups
	for x in rej_susie_meta:
		group = set([orig2new(j) for j in x['group']])
		# We don't consider contiguous groups
		if len(group) - 1 != max(group) - min(group):
			cand_groups.append(
				CandidateGroup(
					group=group, 
					pep=x['pep'],
					data=dict(snps=x['snps']),
				)
			)

	# Run for both FDR
	print(f"Preprocessing finished, running BLiP at time={elapsed(time0)}.")
	output = {}
	for error in ['fdr']:
		input_cand_groups = [copy.deepcopy(x) for x in cand_groups]
		rejections_blip = pyblip.blip.BLiP(
			cand_groups=input_cand_groups,
			weight_fn='inverse_size',
			q=q,
			error=error,
			deterministic=True,
		)
		print(f"Finished with blip rejections for error={error} at time={elapsed(time0)}.")

		# Carefully re-index back to original alphas index
		for rej in rejections_blip:
			rej.data.pop("blip-group")
			rej.data['group'] = list(rej.group) # makes it serializable
			rej.data['pep'] = rej.pep
			# Calculate the set of SNPs. This involves careful indexing.
			# 1. Find the SNPs for indices which are part of 'rel_inds', 
			# as signalled by the fact that the index is less than sum(rel_inds).
			rel_in_group = [j for j in rej.group if j < nrel]
			snps = alphas_index[rel_inds][rel_in_group].tolist()
			# 2. Find SNPs for indices that are not part of `rel_inds', 
			# signalled by the fact that the index is >= sum(rel_inds).
			not_rel_in_group = [int(j - nrel) for j in rej.group if j >= nrel]
			snps = snps + alphas_index[not_rel_in_group].tolist()
			# Do some checks to make sure we got this right
			if len(snps) != len(set(snps)):
				print(f"snps {snps} contains duplicates, which indicates an indexing error.")

			# For groups rejected by susie, we can check that we got the same result
			if 'snps' in rej.data:
				if set(snps) != set(rej.data['snps']):
					raise ValueError(f"snps do not line up {snps} vs. {rej.data['snps']}")
			else:
				rej.data['snps'] = alphas_index[rel_inds][list(rej.group)].tolist()

		# Add to output
		output[error] = [rej.data for rej in rejections_blip]

	# Return
	print(f"Res-adj power for blip: {sum([1/len(x['group']) for x in output['fdr']])}")
	return output

def susie_indiv_results(single_peps, q):
	"""
	Lines up PEPs and rejects as many as possible
	while controlling Bayesian FDR
	"""
	inds = np.argsort(single_peps)
	fdrs = np.cumsum(single_peps[inds]) / np.arange(1, single_peps.shape[0] + 1)
	if len(fdrs) == 0:
		return []
	if fdrs[0] > q:
		return []
	else:
		num_selections = np.where(fdrs <= q)[0].max()
		return [[x] for x in inds[:num_selections]]

def susie_results(alphas, alphas_index, q, max_size, time0):

	# Loop through loci (chromosome x start)
	alphas['rangeindex'] = np.arange(alphas.shape[0]) 
	chromes = alphas['CHR'].unique()
	susie_rejections = []
	for chrome in chromes:
		chrome_sub = alphas.loc[(alphas['CHR'] == chrome)]
		starts = chrome_sub['REGION_START'].unique().tolist()
		for start in tqdm(starts):
			subset = alphas.loc[
					(alphas['CHR'] == chrome) &
					(alphas['REGION_START'] == start)
			]
			rangeind = subset['rangeindex']
			for col in [f'ALPHA{i}' for i in range(1, 11)]:
				alphajs = subset[col].values
				sortinds = np.argsort(-1*alphajs)
				cumsums = np.cumsum(alphajs[sortinds[0:max_size]])
				if cumsums[-1] > 1 - q:
					rej_length = np.where(cumsums > 1 - q)[0].min()
					rej_inds = rangeind[sortinds[0:(rej_length+1)]]
					susie_rejections.append({
						"group":rej_inds.tolist(),
						"pep":1 - cumsums[rej_length],
						"snps":alphas_index[rej_inds.tolist()].tolist()
					})
		print(f"For susie, finished chrome={chrome}, time={elapsed(time0)}.")


	return susie_rejections


############ GLOBALS ################
OUTCOMES = [
	('disease_CARDIOVASCULAR', False),
	('biochemistry_HDLcholesterol', True),
	('biochemistry_LDLdirect', False),
	('body_HEIGHTz', True),
]

MAX_SIZE = 25
DEBUG = False
FUNCT = False
q = 0.05
thresh = 0.25
rel_thresh = 0.99

def main():

	time0 = time.time()

	# Loop through outcomes
	for (outcome, by_chromosome) in OUTCOMES:
		print("==================================================================")
		print("==================================================================")
		print(f"Now analyzing outcome = {outcome}")
		print("==================================================================")
		print("==================================================================")

		# Load the set of snps which we use for MCMC
		all_snps = set()
		range_end = 22 + 1
		if DEBUG:
			range_end = 2

		for chrome in list(range(1, range_end)):
			XTY_snps = set(pd.read_csv(
				f"main_cache/XTY/chrome{chrome}_trait{outcome}.csv",
				index_col=0,
				usecols=[0]
			).index.tolist())
			all_snps = all_snps.union(XTY_snps)
			print(f"For outcome={outcome}, after chrome={chrome}, num_snps={len(all_snps)} at {elapsed(time0)}.")

		# Load data
		alphas, alphas_np = load_alphas(outcome, funct=FUNCT, time0=time0)
		
		# To speed up debugging
		if DEBUG:
			alphas = alphas.iloc[0:10000000]
			alphas_np = alphas_np[0:10000000]

		# Possibly chunk into chromosomes
		if by_chromosome:
			del alphas_np
			alphas_list = []
			alphas_np_list = []
			chromes = alphas['CHR'].unique()
			for chrome in chromes:
				chrome_alphas = alphas.loc[alphas['CHR'] == chrome]
				alphas_list.append(chrome_alphas)
				alphas_np_list.append(chrome_alphas[NUM_COLS].values)
		else:
			alphas_list = [alphas]
			alphas_np_list = [alphas_np]

		# For each chunk, run different analyses
		all_susie_meta = []
		all_blip_meta_fdr = []
		all_susie_indv_meta = []
		all_df = pd.DataFrame(columns=['size', 'method'])
		for alphas, alphas_np in zip(alphas_list, alphas_np_list):

			# SNP names that line up with the MCMC analysis
			alphas['snp'] = alphas['CHR'].astype(str).copy()
			alphas['snp'] += "." + alphas['BP'].astype(str)
			alphas['snp'] += "." + alphas['A2'].astype(str)
			alphas['snp'] += "." + alphas['A1'].astype(str)

			# Only consider SNPs which we run MCMC on for fair comparison
			flags = alphas['snp'].isin(all_snps)
			alphas = alphas.loc[flags]
			alphas_np = alphas_np[flags]

			# PEPs for each feature
			single_peps = 1 - alphas['PIP'].values

			# Run different analyses, collect metadata + indices
			# 1. susie by itself
			rej_susie_meta = susie_results(
				alphas=alphas,
				alphas_index=alphas['snp'].values,
				q=q,
				max_size=MAX_SIZE,
				time0=time0
			)
			all_susie_meta.extend(rej_susie_meta)
			rej_susie = [list(x['group']) for x in rej_susie_meta]
			# 2. BLiP + susie
			rej_blip_meta = blip_results(
				alphas=alphas_np,
				alphas_index=alphas['snp'].values,
				single_peps=single_peps,
				max_size=MAX_SIZE,
				rej_susie_meta=rej_susie_meta,
				q=q,
				thresh=thresh,
				rel_thresh=rel_thresh,
				time0=time0
			)
			all_blip_meta_fdr.extend(rej_blip_meta['fdr'])
			rej_blip_fdr = [x['group'] for x in rej_blip_meta['fdr']]
			# 3. susie only with individual PEPs
			rej_susie_indv = susie_indiv_results(
				single_peps=single_peps, q=q
			)
			snpvals = alphas['snp'].values
			rej_susie_indv_meta = [
				dict(
					snps=snpvals[x].tolist(),
					pep=single_peps[x][0]
				) for x in rej_susie_indv
			]
			all_susie_indv_meta.extend(rej_susie_indv_meta)

			# Format results
			df_list = []
			for rej, method in zip(
				[rej_blip_fdr, rej_susie, rej_susie_indv],
				['susie_blip', 'susie', 'susie-indiv-only']
			):

				rej_size = np.array([len(x) for x in rej])
				df = pd.DataFrame(columns=['size', 'method'])
				df['size'] = rej_size
				df['method'] = method
				df_list.append(df)
			df = pd.concat(df_list, axis='index')

			# Join with other chunk (chromosome) data
			all_df = pd.concat([all_df, df], axis='index')
			all_df.to_csv(f"output/polyfun_wrapper/groupsizes/{outcome}.csv")

			# Save output metadata
			output_dict = {
				"susie_blip":all_blip_meta_fdr,
				"susie":all_susie_meta,
				"susie-indiv-only":all_susie_indv_meta,
			}
			with open(f"output/polyfun_wrapper/rejections/{outcome}.json", 'w') as thefile:
				thefile.write(json.dumps(output_dict))



if __name__ == '__main__':
	main()