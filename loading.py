import sys
import os

import numpy as np
import pandas as pd
import scipy.sparse


def load_ld(chrome, start, end, return_df=True):
	# Read LD data, start with metadata
	loci = f"chr{chrome}_{start}_{end}"
	meta = pd.read_table(f"data/ld/{loci}.gz", sep="\s+")
	meta.index = meta['chromosome'].astype(str) + "." + meta['position'].astype(str)
	meta.index = meta.index + "." + meta['allele2'].astype(str)

	# Actual correlation matrix
	ld = scipy.sparse.load_npz(f"data/ld/{loci}.npz").toarray()
	ld += ld.T

	# Turn into df
	if return_df:
		ld = pd.DataFrame(ld, index=meta.index, columns=meta.index)

	return (ld, meta)

def load_sumstats(outcome, bolt_lmm_inf=False, all_cols=False):

	# Determine which cols to read
	if all_cols:
		names = None
	else:
		names = ["CHR", "BP", "ALLELE0", "CHISQ_BOLT_LMM", "BETA", "SE"]

	# Read data
	filename = f"bolt_337K_unrelStringentBrit_MAF0.001_v3.{outcome}"
	sumstats = pd.read_table(
		f"data/polyfun_results/{filename}.bgen.stats.gz",
		sep="\s+",
		names=names
	)
	sumstats.index = sumstats['CHR'].astype(str) + "." + sumstats['BP'].astype(str)
	sumstats.index = sumstats.index + "." + sumstats["ALLELE0"].astype(str)

	# Filter out errors where chisq < 0
	sumstats = sumstats.loc[
		sumstats['CHISQ_BOLT_LMM'] > 0
	]

	# Compute Z-statistics
	if bolt_lmm_inf:
		sumstats['Z'] = sumstats['BETA'] / sumstats['SE']
	else:
		sumstats['Z'] = np.sqrt(
			sumstats['CHISQ_BOLT_LMM']
		) * np.sign(sumstats['BETA'])

	# Return pd.series
	return sumstats

def intersect_ld_sumstats(sumstats, ld):

	# Find intersection of LD/sumstat snps
	ss_snps = set(sumstats.index.tolist())
	ld_snps = set(ld.index.tolist())
	int_snps = ld_snps.intersection(ss_snps)

	# Subset summary statistics + ld info
	sstats = sumstats.loc[int_snps]
	ldmat = ld.loc[int_snps, int_snps]

	return sstats, ldmat


def run_susie(sumstats, ld, n, L=10):
	"""
	Runs susie.
	Code adapted from https://github.com/omerwe/polyfun/blob/master/finemapper.py
	"""

	# Get XTX and XTy
	sstats, ldmat = intersect_ld_sumstats(sumstats, ld)

	#load SuSiE R package
	import rpy2
	import rpy2.robjects.numpy2ri as numpy2ri
	import rpy2.robjects as ro
	ro.conversion.py2ri = numpy2ri
	numpy2ri.activate()
	from rpy2.robjects.packages import importr
	susieR = importr('susieR')
	R_null = ro.rinterface.NULL
	#self.RNULLType = rpy2.rinterface.RNULLType

	#rpy2 bug fix
	import rpy2.robjects.numpy2ri as numpy2ri
	reload(numpy2ri)
	numpy2ri.activate()

	# Run susie
	p = ldmat.shape[0]
	Zstats = sstats['Z'].values.copy()
	susie_obj = susieR.susie_suff_stat(
		bhat=Zstats.reshape((-1,1)),
		shat=np.ones((p,1)),
		R=ldmat.values,
		n=n,
		L=L,
		scaled_prior_variance=0.0001,
		estimate_prior_variance=True,
		residual_variance=R_null,
		estimate_residual_variance=True,
		verbose=True,
		prior_weights=R_null
	)
	return susie_obj
	