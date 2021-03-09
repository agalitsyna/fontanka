
# import standard python libraries
import glob
import multiprocess
import re

# data analysis
import numpy as np
import scipy
import pandas as pd
from scipy import signal

# plotting libraries
import matplotlib.pyplot as plt
import seaborn as sns
import proplot

# ignore warnings
import warnings
warnings.filterwarnings('ignore')

# import libraries for Hi-C data analysis
import cooler
import bioframe

import cooltools
import cooltools.expected
from cooltools import snipping
import cooltools.lib.plotting


# Additional imports:
from cooltools.lib._query import CSRSelector
from cooltools.lib import peaks, numutils
from cooltools.insulation import get_n_pixels

SCHARR = np.array([[-3 - 3j, 0 - 10j, +3 - 3j],
                   [-10 + 0j, 0 + 0j, +10 + 0j],
                   [-3 + 3j, 0 + 10j, +3 + 3j]])  # Gx + j*Gy


### Insulation-based fountain score (outdated):
def _find_insulating_boundaries_dense(
        clr,
        window_bp=100000,
        balance="weight",
        min_dist_bad_bin=2,
        ignore_diags=None,
        chromosomes=None,
        kernel=None,
        mode='min',
        with_shar_filter=True
):
    """Calculate the diamond insulation scores and call insulating boundaries.
    Parameters
    ----------
    c : cooler.Cooler
        A cooler with balanced Hi-C data.
    window_bp : int
        The size of the sliding diamond window used to calculate the insulation
        score.
    min_dist_bad_bin : int
        The minimal allowed distance to a bad bin. Do not calculate insulation
        scores for bins having a bad bin closer than this distance.
    ignore_diags : int
        The number of diagonals to ignore. If None, equals the number of
        diagonals ignored during IC balancing.
    kernel : numpy.array
        A dense matrix of kernel
    mode : string, 'min' (boundaries) or 'max' (dense regions)
        Find the maximum or minimum peaks

    Returns
    -------
    ins_table : pandas.DataFrame
        A table containing the insulation scores of the genomic bins and
        the insulating boundary strengths.
    """
    if chromosomes is None:
        chromosomes = clr.chromnames

    bin_size = clr.info["bin-size"]
    ignore_diags = (
        ignore_diags
        if ignore_diags is not None
        else clr._load_attrs(clr.root.rstrip("/") + "/bins/weight")["ignore_diags"]
    )
    window_bins = window_bp // bin_size

    if window_bp % bin_size != 0:
        raise Exception(
            "The window size ({}) has to be a multiple of the bin size {}".format(
                window_bp, bin_size
            )
        )

    ins_chrom_tables = []
    for chrom in chromosomes:
        ins_chrom = clr.bins().fetch(chrom)[["chrom", "start", "end"]]
        is_bad_bin = np.isnan(clr.bins().fetch(chrom)["weight"].values)

        m = clr.matrix(balance=balance).fetch(chrom)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            ins_run = _insul_diamond_dense(m, window_bins, ignore_diags, kernel=kernel,
                                           with_shar_filter=with_shar_filter)
            ins_track = ins_run[0]
            pixels_track = ins_run[1]
            scharr_track = ins_run[3]
            ins_track = np.log2(ins_track)

        bad_bin_neighbor = np.zeros_like(is_bad_bin)
        for i in range(0, min_dist_bad_bin):
            if i == 0:
                bad_bin_neighbor = bad_bin_neighbor | is_bad_bin
            else:
                bad_bin_neighbor = bad_bin_neighbor | np.r_[[True] * i, is_bad_bin[:-i]]
                bad_bin_neighbor = bad_bin_neighbor | np.r_[is_bad_bin[i:], [True] * i]

        ins_track[bad_bin_neighbor] = np.nan
        ins_chrom["bad_bin_masked"] = bad_bin_neighbor

        ins_track[~np.isfinite(ins_track)] = np.nan

        ins_chrom["log2_insulation_score_{}".format(window_bp)] = ins_track

        if mode == 'min':
            poss, proms = peaks.find_peak_prominence(-ins_track)
        elif mode == 'max':
            poss, proms = peaks.find_peak_prominence(ins_track)
        else:
            raise Exception(f'Mode: {mode} is not implemented.')

        ins_prom_track = np.zeros_like(ins_track) * np.nan
        ins_prom_track[poss] = proms
        ins_chrom["boundary_strength_{}".format(window_bp)] = ins_prom_track
        ins_chrom["n_valid"] = pixels_track
        ins_chrom["scharr_value"] = scharr_track

        ins_chrom_tables.append(ins_chrom)

    ins_table = pd.concat(ins_chrom_tables)
    return ins_table

def _insul_diamond_dense(mat,
                         window=10,
                         ignore_diags=2,
                         norm_by_median=True,
                         kernel=None,
                         with_shar_filter=True):
    """
    Calculates the insulation score of a Hi-C interaction matrix.
    Parameters
    ----------
    mat : numpy.array
        A dense square matrix of Hi-C interaction frequencies.
        May contain nans, e.g. in rows/columns excluded from the analysis.
    window : int
        The width of the window to calculate the insulation score.
    ignore_diags : int
        If > 0, the interactions at separations < `ignore_diags` are ignored
        when calculating the insulation score. Typically, a few first diagonals
        of the Hi-C map should be ignored due to contamination with Hi-C
        artifacts.
    norm_by_median : bool
        If True, normalize the insulation score by its NaN-median.
    kernel : numpy.array
        A kernel matrix used for convolution with insulation window.
    """
    if ignore_diags:
        mat = mat.copy()
        for i in range(-ignore_diags + 1, ignore_diags):
            numutils.set_diag(mat, np.nan, i)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        N = mat.shape[0]
        sum_balanced = np.nan * np.ones(N)
        n_pixels = np.zeros(N)
        scharr_value = np.nan * np.ones(N)
        for i in range(0, N):
            lo = max(0, i + 1 - window)
            hi = min(i + window, N)
            # nanmean of interactions to reduce the effect of bad bins
            snip = mat[lo: i + 1, i:hi]
            if not kernel is None:
                lo_kern = max(0, window - (i + 1))
                hi_kern = min(window, N - i)
                snip = np.multiply(kernel[lo_kern:, :hi_kern], snip)
            sum_balanced[i] = np.nanmean(snip)
            n_pixels[i] = np.sum(np.isfinite(snip))

            if with_shar_filter:
                grad = signal.convolve2d(snip, SCHARR, boundary='symm', mode='same')
                if not kernel is None:
                    grad = np.multiply(kernel[lo_kern:, :hi_kern], grad)
                scharr_value[i] = np.nanmax(np.absolute(grad))

        score = sum_balanced / n_pixels

        if norm_by_median:
            score /= np.nanmedian(score)

    return score, n_pixels, sum_balanced, scharr_value

### Tools for snip-based insulation
def generate_ObsExpSnips(clr, windows, regions, nthreads=10):

    # Calculate expected interactions for chromosome arms
    with multiprocess.Pool(nthreads) as pool:
        expected = cooltools.expected.diagsum(
            clr,
            regions=regions,
            transforms={
                'balanced': lambda p: p['count'] * p['weight1'] * p['weight2']
            },
            map=pool.map
        )
    expected = pd.concat(
        [expected, bioframe.region.parse_regions(expected.region)],
        axis=1).drop('name', axis=1)
    expected['balanced.avg'] = expected['balanced.sum'] / expected['n_valid']

    snipper = cooltools.snipping.ObsExpSnipper(clr, expected)
    with multiprocess.Pool(nthreads) as pool:
        stack = cooltools.snipping.pileup(
                windows,
                snipper.select,
                snipper.snip,
                map=pool.map
                )
    return stack

def save_snips(stack, outfile):
    with open(outfile, 'wb') as f:
        np.save(f, stack)

def read_snips(infile):
    with open(infile, 'rb') as f:
        stack = np.load(f)
    return stack

def generate_fountain_score(stack, kernel):
    n_snips = stack.shape[2]
    #snip_size = stack.shape[0]
    track_fs = np.ones( n_snips )*np.nan
    for i in range(n_snips):
        snip = stack[:, :, i]
        fs = np.multiply(kernel, snip)
        track_fs[i] = np.nanmean(fs)
    return track_fs

def generate_scharr_score(stack, kernel=None):
    n_snips = stack.shape[2]
    #snip_size = stack.shape[0]
    track_scharr = np.ones( n_snips )*np.nan
    for i in range(n_snips):
        snip = stack[:, :, i].copy()
        if not kernel is None:
            snip = snip * kernel
        grad = signal.convolve2d(snip, SCHARR, boundary='symm', mode='same')
        scharr = np.nanmax( np.absolute(grad) )
        track_scharr[i] = scharr
    return track_scharr

def get_peaks_prominence(track):
    """ Very simple wrapper to get the prominence of peaks (no filtering)"""
    poss, proms = peaks.find_peak_prominence(track)
    ins_prom_track = np.zeros_like(track) * np.nan
    ins_prom_track[poss] = proms
    return ins_prom_track

### Tools for Fountains and eigenvectors analysis of snips:

def inverted_trinagle_mask(alpha, h, fill_pos=1, fill_neg=0):
    """
    Return square mask with a single fountain (from lower left to upper right).
    """
    # No obtuse angles are allowed
    mtx = np.zeros([h,h])+fill_neg
    for x in range(h):
        bgn = np.trunc((x+1)*np.tan((np.pi/2-alpha)/2)).astype(int)
        end = np.ceil((x+0.5)/np.tan((np.pi/2-alpha)/2)).astype(int)
        bgn = bgn if bgn>0 else 0
        end = end if end<h else h
        mtx[x, bgn:end] = fill_pos
    # if norm_sum:
    #     return np.flip( mtx/np.sum(mtx) , axis=1).T
    return np.flip( mtx , axis=1).T

def double_triangle_mask(alpha, size, fill_pos=1, fill_neg=0):
    """
    Return square mask with double fountains. The h is the half-size. The resulting size will be 2*h + 1 .
    """
    mask_triangle = inverted_trinagle_mask(alpha, size, fill_pos=fill_pos, fill_neg=fill_neg)
    mask_double_triangle = np.vstack([
        np.hstack([np.ones([size, size]) * (fill_neg), mask_triangle]),
        np.hstack([np.flip(mask_triangle), np.ones([size, size]) * (fill_neg)])
    ])
    mask = (fill_neg) * np.ones([2 * size + 1, 2 * size + 1])
    mask[np.triu_indices(2 * size + 1, 1)] = mask_double_triangle[np.triu_indices(2 * size, 0)]
    mask[np.tril_indices(2 * size + 1, -1)] = mask_double_triangle[np.tril_indices(2 * size, 0)]

    return mask

def get_eig(snip, n_eigs):
    snip = snip.copy()
    snip[np.isnan(snip)] = 0
    eigvecs, eigvals = numutils.get_eig(snip, n_eigs, mask_zero_rows=True)
    eigvecs /= np.sqrt(np.nansum(eigvecs ** 2, axis=1))[:, None]
    eigvecs *= np.sqrt(np.abs(eigvals))[:, None]
    return eigvecs

def flip_eigs(eigs, vect):
    """ Flipping the sign accoring to the vector sign. """
    return (eigs.T*np.sign(vect)).T


def get_components(stack, n_eigs):
    n_snips = stack.shape[2]
    snip_size = stack.shape[0]
    stack_eigvects = np.ones((n_eigs, snip_size, n_snips)) * np.nan
    for i in range(n_snips):
        snip = stack[:, :, i]
        snip[np.isnan(snip)] = 0
        if np.sum(snip) == 0:
            continue
        try:
            eigvecs, eigvals = numutils.get_eig(snip, n_eigs + 1, mask_zero_rows=True)
            eigvecs /= np.sqrt(np.nansum(eigvecs ** 2, axis=1))[:, None]
            eigvecs *= np.sqrt(np.abs(eigvals))[:, None]

            stack_eigvects[:, :, i] = eigvecs[1:, :]
        except Exception:
            pass
    return stack_eigvects

def flip_sign(eigvects, n_eigs, reference_eigvects):
    n_snips   = eigvects.shape[2]
    snip_size = eigvects.shape[1]
    stack_eigvects = np.ones((n_eigs, snip_size, n_snips)) * np.nan
    for i in range(n_snips):
        eigs = eigvects[:, :, i]
        eigs[np.isnan(eigs)] = 0
        if np.sum(eigs) == 0:
            continue
        try:
            vect = np.array( [scipy.stats.pearsonr(eigs[i, :], reference_eigvects[i, :])[0] for i in range(n_eigs)] )
            stack_eigvects[:, :, i] = flip_eigs(eigs, vect)
        except Exception:
            pass
    return stack_eigvects








