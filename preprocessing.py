import sys
import os
import warnings

import time
import copy
import random
import numpy as np
import scipy as sp
import networkx as nx
import networkx.algorithms.mis as mis

import wget
import requests

def url_exists(url):
    if requests.head(url).status_code == 200:
        return True
    else:
        return False

def elapsed(time0):
    return np.around(time.time() - time0, 2)

def min_eigval(cov):
    """
    eigsh is faster for super high dimensions,
    but less accurate and reliable
    and slower in low dimensions?
    """
    #     try:
    #         return eigsh(cov, 1, which='SA')[0][0]
    #     except:
    #        return np.linalg.eigh(cov)[0].min()
    return np.linalg.eigh(cov)[0].min()

def shift_until_PSD(
    cov, 
    tol, 
    n_iter=8,
    time0=None,
    init_gamma=None,
    conv_tol=1 / (2**10),
):
    p = cov.shape[0]
    if p < 7500:
        mineig = min_eigval(cov)
        if mineig < tol:
            gamma = (tol - mineig) / (1 - mineig)
        else:
            gamma = 0
        return cov * (1-gamma) + gamma * np.eye(p)
    else:
        time0 = time.time() if time0 is None else time0
        ugamma = 0.2 # min gamma controlling eig bound, if > 0.3 this is really bad
        lgamma = 0 # max gamma violating eig bound
        for j in range(n_iter):
            if init_gamma is not None and j == 0:
                gamma = init_gamma
            else:
                gamma = (ugamma + lgamma) / 2
            try:
                print(f"Trying cholesky with gamma={gamma} for p={p} at t={elapsed(time0)}")
                sys.stdout.flush()
                np.linalg.cholesky(cov * (1 - gamma) + (gamma - tol) * np.eye(p))
                ugamma = gamma
            except np.linalg.LinAlgError:
                lgamma = gamma
            print(f"After it={j}, ugamma={ugamma}, lgamma={lgamma}")
            if ugamma - lgamma < conv_tol:
                break
        return cov * (1 - ugamma) + ugamma * np.eye(p)

def force_PSD_mis(
    ld, 
    time0=None, 
    max_corr=0.999,
    max_cset=1,
    tol=1e-5
):
    time0 = time.time() if time0 is None else time0
    # For reproducability (max ind set is random alg)
    np.random.seed(12345)
    random.seed(12345)
    # Find maximal independent set
    p = ld.shape[0]
    flags = np.abs(ld) >= max_corr
    for j in range(flags.shape[0]):
        flags[j, j] = False
    G = nx.Graph(flags)
    mis_inds = sorted(list(mis.maximal_independent_set(G)))
    # possibly shift until psd
    M = shift_until_PSD(
        ld[mis_inds][:, mis_inds],
        tol=tol, 
        init_gamma=0.025,
        time0=time0
    )
    for kind, k in enumerate(mis_inds):
        ld[k][mis_inds] = M[kind]
    
    # adjust the other columns of ld
    to_adjust = set(range(p)) - set(mis_inds)
    adjusted_inds = copy.copy(mis_inds)
    l_nadj = sorted(list(to_adjust))
    mu_mat = np.zeros((len(to_adjust), p))
    
    for jind, j in enumerate(l_nadj):
        # variables most highly corr with j
        maxcorrinds = np.argsort(-1*np.abs(ld[j]))
        cset = []
        for k in maxcorrinds:
            if k in adjusted_inds:
                if len(cset) == 0:
                    cset.append(k)
                elif np.max(np.abs(ld[k][cset])) < 0.99:
                    cset.append(k)
            if len(cset) >= max_cset:
                break
    
        # form cov matrix and submatrices
        cset_pj = copy.copy(cset)
        cset_pj.insert(0, j)
        cov = ld[cset_pj][:, cset_pj]
        cov = shift_until_PSD(cov, tol=tol) # make sure PSD
        covnegj = cov[1:][:, 1:]
        covjnegj = cov[0][1:]
        dot_j_negj = sp.linalg.solve(
            a=covnegj,
            b=covjnegj,
            assume_a='pos'
        )
        # set of linear transformations
        mu_mat[jind][cset_pj[1:]] = dot_j_negj
        # print progress
        if jind % 1000 == 0:
            print(f"Finished with jind={jind} at {elapsed(time0)}")
    
    # Calculate full covariance matrix
    print(f"Calculating full matrix at {elapsed(time0)}")
    mu_ld = np.dot(mu_mat, ld) # k x p
    nadj2nadj = np.dot(mu_ld, mu_mat.T)
    nadj2nadj += np.diag(1 - np.diag(nadj2nadj))
    l_adj = sorted(list(adjusted_inds))
    output = np.concatenate(
        [np.concatenate(
            [ld[l_adj][:, l_adj], mu_ld[:, l_adj]], axis=0
        ),
         np.concatenate(
             [mu_ld[:, l_adj].T, nadj2nadj], axis=0
        )],
        axis=1
    )
    # re-order output
    print(f"Reordering ld matrix at {elapsed(time0)}")
    rev_inds = np.zeros(p).astype(int)
    for i, j in enumerate(l_adj + l_nadj):
        rev_inds[j] = i
    r = output[rev_inds][:, rev_inds]
    return shift_until_PSD(r, tol=tol, init_gamma=0.0025, time0=time0)

def download_ld_data(chrome, start):

    # Set directory to data/ld/
    file_directory = os.path.dirname(os.path.abspath(__file__))
    os.chdir(f"{file_directory}/data/ld/")

    # find urls and wget
    end = int(start+3000000)
    baseurl = f"https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/chr{chrome}"
    for postfix in ['gz', 'npz']:
        url = f"{baseurl}_{int(start)}_{end}.{postfix}"
        fname = f"chr{chrome}_{start}_{end}.{postfix}"
        if not os.path.exists(fname):
            wget.download(url)

    # Reset working directory
    os.chdir(file_directory)
    return 0

def modified_HESS(
    ld,
    sumstats,
    n,
    max_corr=0.95,
    hess_iter=1,
    prop_keep=0.005,
    max_hess_size=10000
):
    hg2s = []
    p = sumstats.shape[0]
    sumstats_cutoff = np.quantile(np.abs(sumstats), 1-prop_keep)
    inds = np.where(np.abs(sumstats) > sumstats_cutoff)[0]
    if len(inds) == 0:
        inds = list(range(p))
    # preprocessing for mis
    ld_sub = ld[inds][:, inds]
    flags = np.abs(ld_sub) >= max_corr
    for j in range(flags.shape[0]):
        flags[j, j] = False
    G = nx.Graph(flags)
    # Loop through and create estimates
    for _ in range(hess_iter):
        mis_inds = sorted(list(mis.maximal_independent_set(G)))
        m = len(mis_inds)
        if m > max_hess_size:
            mis_inds = np.random.choice(m, max_hess_size, replace=False)
        sumstats_sub = sumstats[inds][mis_inds]
        alphaRinv = sp.linalg.solve(
            ld_sub[mis_inds][:, mis_inds],
            sumstats_sub,
            assume_a='pos'
        )
        hg2s.append(np.dot(sumstats_sub, alphaRinv) - m / n)

    return np.mean(np.array(hg2s))