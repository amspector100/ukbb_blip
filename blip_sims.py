import time
import json
import datetime
import copy
import numpy as np
import pandas as pd
import loading
import warnings
import scipy.linalg
import parser
import networkx as nx
import networkx.algorithms.mis as mis

import os
import sys
file_directory = os.path.dirname(os.path.abspath(__file__))
parent_directory = os.path.split(file_directory)[0]
sys.path.insert(0, parent_directory + "/pyblip/")
import pyblip
print(f"pyblip version is {pyblip.__version__}")

# Local modules
import utilities
import preprocessing
from preprocessing import elapsed, shift_until_PSD, min_eigval

COLUMNS = ['power', 'fdr', 'method', 'hg2', 'num_causal', 'seed']

def load_processed_ld(chrome, start, time0):
	# Check if processed ld is cached
	end = int(start+3000000)
	processed_fname = f'sims/processed_ld/chr{chrome}_{start}_{end}.npy'
	if os.path.exists(processed_fname):
		return np.load(processed_fname)
	else:
		# Make sure LD data is downloaded
		preprocessing.download_ld_data(
			chrome=chrome, start=start
		)
		# Load ld
		ld_file = f"data/ld/chr{chrome}_{start}_{end}.npz"
		ld = scipy.sparse.load_npz(ld_file).toarray().astype(np.float32)
		ld += ld.T
		ld = preprocessing.force_PSD_mis(
			ld, 
			time0=time0,
			tol=1e-5
		).astype(np.float32)
		# Cache
		np.save(processed_fname, ld)
		return ld

def sample_XTY_noise(
	L,
	reps,
	max_corr=0.99,
	tol=1e-3,
	time0=None,
	max_cset=1,
	cache=True,
):
	"""
	L is the cholesky decomposition of the LD matrix.
	(This function assumes the LD matrix has already
	been preprocessed using the MIS technique from 
	from Weisbrod et al. (2019), so it is PSD.)
	"""
	# Initialize samples from N(0, ld)
	p = L.shape[0]
	samples = np.random.randn(p, reps)
	samples = np.dot(
		L, samples
	)

	return samples

def create_XTY(n, noisei, ld, hg2, num_causal):
	p = ld.shape[0]
	# Create sparse coefficients
	beta = np.zeros(p)
	causal_inds = np.random.choice(
		np.arange(p),
		num_causal,
		replace=False,
	)
	beta[causal_inds] = np.random.randn(num_causal)

	# Adjust beta so that hg2 = 2
	scale = np.dot(beta, np.dot(ld, beta))
	beta = np.sqrt(hg2) * beta / np.sqrt(scale)
	
	# Create XTY where var(Y) = 1
	XTY = np.dot(ld, beta) + np.sqrt((1 - hg2)/n) * noisei
	XTY = n * XTY
	return XTY, beta

# For running SuSiE
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
from rpy2.rinterface_lib.embedded import RRuntimeError
def run_susie_suff_stats(
	XTY,
	ld,
	n,
	L,
	q,
	hg2_est,
	**kwargs
):
	# This code adapted from the polyfun package
	# load SuSiE R package
	ro.conversion.py2ri = numpy2ri
	numpy2ri.activate()
	from rpy2.robjects.packages import importr
	susieR = importr('susieR')
	R_null = ro.rinterface.NULL
	# Run susie
	prior_variance = hg2_est / L
	susie_obj = susieR.susie_suff_stat(
		Xty=XTY.reshape(-1, 1),
		XtX=n * ld,
		n=n,
		L=L,
		yty=n,
		scaled_prior_variance=prior_variance,
		estimate_prior_variance=False,
		#residual_variance=1-hg2_est,
		estimate_residual_variance=False,
		**kwargs
	)
	# Extract output
	alphas = np.array(susie_obj.rx2('alpha'))
	susie_sets = susie_obj.rx2('sets')[0]
	try:
		susie_sets = [
			np.array(s)-1 for s in susie_sets
		]
	except TypeError:
		susie_sets = []

	return alphas, susie_sets

def run_susie_blip(
	seed,
	n,
	noise,
	hg2,
	num_causal,
	ld,
	args,
	time0,
):
	# For logging
	dgp_args = [hg2, num_causal, seed]
	print_id = f"seed={seed}, hg2={hg2}, num_causal={num_causal}"

	# Step 1: sample XTY and create true coeffs
	XTY, beta = create_XTY(
		n=n,
		noisei=noise[:, seed],
		ld=ld, 
		hg2=hg2,
		num_causal=num_causal,
	)

	# Step 2: estimate hg2 using modified HESS
	print(f"Estimating hg2 via HESS for {print_id} at {elapsed(time0)}.")
	hg2_est = preprocessing.modified_HESS(
		ld=ld,
		sumstats=XTY / n,
		n=n,
		max_corr=args.get("max_corr_hess", [0.95])[0],
		hess_iter=args.get("hess_iter", [10])[0],
	)

	# Step 3: fit SuSiE
	print(f"Starting SuSiE with hg2_est={hg2_est} for {print_id} at {elapsed(time0)}.")
	hg2_est = max(hg2_est, 1e-6)
	# Note alphas is L x p
	q=args.get("q", [0.05])[0]
	alphas, susie_sets = run_susie_suff_stats(
		XTY=XTY,
		ld=ld,
		n=n,
		L=args.get("L", [10])[0],
		hg2_est=hg2_est,
		q=q,
		verbose=args.get("verbose", [False])[0],
	)
	print(f"Finished fitting SuSiE for {print_id} at {elapsed(time0)}.")

	# Step 4: run BLiP after prefiltering SNPs 
	# to exclude SNPs with marginal PIP < 0.01
	prefilter_threshold = args.get("prefilter_threshold", [0.01])[0]
	marg_pips = 1 - np.exp(np.sum(np.log(np.maximum(1-alphas, 1e-10)), axis=0))
	rel_inds = sorted(np.where(marg_pips > prefilter_threshold)[0])
	if len(rel_inds) == 0:
		detections = []
	else:
		# Create cand groups
		max_pep = args.get("max_pep", [0.25])[0]
		cand_groups = pyblip.create_groups.susie_groups(
			alphas=alphas[:, rel_inds],
			X=None,
			max_pep=max_pep,
			q=q,
			max_size=args.get("max_size", [25])[0],
		)
		# BLiP detections
		detections = pyblip.blip.BLiP(
			cand_groups=cand_groups,
			error='fdr',
			q=q,
			max_pep=max_pep,
		)
		# Map detections back to original indices
		for x in detections:
			x.group = set([
				int(rel_inds[j]) for j in x.group
			])

	# Step 5: Compute power and FDR and return
	susie_power, susie_fdr = utilities.rejset_power(
		rej_sets=susie_sets, beta=beta
	)
	blip_power, blip_fdr = utilities.rejset_power(
		rej_sets=[x.group for x in detections], beta=beta
	)
	output = [
		[susie_power, susie_fdr, 'SuSiE'] + dgp_args,
		[blip_power, blip_fdr, 'SuSiE + BLiP'] + dgp_args
	]
	print(output)
	return output

def main(args):

	# for logging
	print("Starting to run main...")
	time0 = time.time()

	# Create output directory
	today = str(datetime.date.today())
	hour = str(datetime.datetime.today().time())
	hour = hour.replace(':','-').split('.')[0]
	output_path = f'sims/results/{today}/{hour}/results.csv'
	args_path = f'sims/results/{today}/{hour}/args.json'
	output_dir = os.path.dirname(output_path)
	output_file = output_dir + "/results.csv"
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	print(f"Output path is {output_path}")

	# globals / args
	args = parser.parse_args(args)
	with open(args_path, "w") as argsfile:
		argsfile.write(json.dumps(args))

	# Extract arguments
	n = args.get("n", [380000])[0]
	tol = args.get("tol", [1e-4])[0]
	reps = args.get("reps", [16])[0]
	rep_start = args.get("rep_start", [0])[0]
	# select loci
	chrome = args.get("chrome", [10])[0]
	start = args.get("start", [134000001])[0]
	end = int(start + 3000000)
	# Following notation in functionally informed mapping
	hg2 = args.get("hg2", [0.05])[0]
	num_causal = args.get('num_causal', [10])[0]

	# Step 1: load ld
	ld = load_processed_ld(chrome=chrome, start=start, time0=time0)
	L = np.linalg.cholesky(ld) # takes ~10 seconds
	print(f'Finished loading LD and performing cholesky decomp at {elapsed(time0)}')

	# Loop through simulation settings
	all_outputs = []
	for hg2 in args.get("hg2", [0.005]):
		for num_causal in args.get('num_causal', [10]):
			# Step 2: sample noise for XTY
			noise = sample_XTY_noise(
				L=L,
				reps=reps,
				max_corr=args.get("max_corr", [0.99])[0],
				tol=tol,
				time0=time0,
			)
			# Step 3: run SuSiE and SuSiE + BLiP
			output = utilities.apply_pool(
				seed=list(range(reps)),
				func=run_susie_blip,
				constant_inputs=dict(
					noise=noise,
					n=n,
					hg2=hg2,
					num_causal=num_causal,
					ld=ld,
					args=args,
					time0=time0,
				),
				num_processes=args.get(
					"num_processes", [1]
				)[0],
			)
			for x in output:
				all_outputs.extend(x)
			# Construct output df and save
			out_df = pd.DataFrame(
				all_outputs, columns=COLUMNS
			)
			out_df.to_csv(output_file)
			print(out_df)
			print(f"Finished with hg2={hg2}, num_causal={num_causal} at {elapsed(time0)}!")





	# #### Step 2: run NPrior + BLiP
	# if args.get("run_blip", [True])[0]:
	# 	all_output = []
	# 	for how_reg in how_regs:
	# 		output = knockpy.utilities.apply_pool(
	# 			func=run_blip,
	# 			num_processes=args.get("num_processes", [6])[0],
	# 			constant_inputs=dict(
	# 				q=args.get("q", [0.05])[0],
	# 				prefix=prefix,
	# 				how_reg=how_reg,
	# 				time0=time0,
	# 				args=args,
	# 				output_dir=output_dir
	# 			),
	# 			i=list(range(rep_start, rep_start+reps))
	# 		)
	# 		for x in output:
	# 			all_output.extend(x)
	# 		out_df = pd.DataFrame(
	# 			all_output,
	# 			columns=[
	# 				"seed",
	# 				"gap",
	# 				"fdp",
	# 				"num_nnulls",
	# 				"n_false_disc",
	# 				"inc_type",
	# 				"how_reg"
	# 			]
	# 		)
	# 		print(out_df)
	# 		out_df.to_csv(output_path)

if __name__ == "__main__":
	main(sys.argv)