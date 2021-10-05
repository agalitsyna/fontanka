# import standard python libraries
import glob
import multiprocess
import re
from tqdm import tqdm

# Packages for data analysis
import numpy as np
import scipy
import pandas as pd
from scipy import signal

# Packages for plotting
import matplotlib.pyplot as plt
import seaborn as sns
import proplot

# Ignore warnings
import warnings

warnings.filterwarnings("ignore")

# Manage logging
from . import get_logger

logger = get_logger(__name__)

# Import libraries for Hi-C data analysis
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

FILTER_SCHARR = np.array(
    [
        [-3 - 3j, 0 - 10j, +3 - 3j],
        [-10 + 0j, 0 + 0j, +10 + 0j],
        [-3 + 3j, 0 + 10j, +3 + 3j],
    ]
)  # Gx + j*Gy

#######################
### Tools for snip-based convolution:
#######################
def generate_ObsExpSnips(clr, windows, regions, expected=None, nthreads=10):
    """
    Generate observed over expected snips for the cool file.
    :param clr: input cooler, already parsed
    :param windows: input windows
    :param regions: regions table, will be replaced with viewframe in the future.
    :param expected: expected table, pandas dataframe
    :param nthreads: number of threads
    :return: stack, numpy array
    """
    # Calculate expected interactions for chromosome arms
    if expected is None:
        logger.info("Generating expected...")
        with multiprocess.Pool(nthreads) as pool:
            expected = cooltools.expected.diagsum(
                clr,
                regions=regions,
                transforms={"balanced": lambda p: p["count"] * p["weight1"] * p["weight2"]},
                map=pool.map,
            )
        expected = pd.concat(
            [expected, bioframe.region.parse_regions(expected.region)], axis=1
        ).drop("name", axis=1)
        expected["balanced.avg"] = expected["balanced.sum"] / expected["n_valid"]


    logger.info("Generating stack of snips...")
    snipper = cooltools.snipping.ObsExpSnipper(clr, expected, regions=regions)
    with multiprocess.Pool(nthreads) as pool:
        stack = cooltools.snipping.pileup(
            windows, snipper.select, snipper.snip, map=pool.map
        )
    return stack


def save_snips(stack, outfile):
    """
    Save stack of snips into the file.
    :param stack: input stack (3D numpy array)
    :param outfile: output file (will be in .npy format)
    :return: None
    """
    logger.info(f"Saving stack into {outfile}...")
    with open(outfile, "wb") as f:
        np.save(f, stack)


def read_snips(infile):
    """
    Read snips from the file.
    :param infile: input file in .npy format
    :return: output stack (3D numpy array)
    """
    logger.info(f"Reading stack from {infile}...")
    with open(infile, "rb") as f:
        stack = np.load(f)
    return stack


def generate_fountain_score(stack, kernel):
    """
    Generate fountain score for the stack.
    :param stack: input stack (3D numpy array)
    :param kernel: kernel (or mask, 2D numpy array) for convolution,
                    size should correspond to the first two dimensions of stack
    :return: track with fountain score (numpy array)
    """
    logger.info(f"Generating fountain score...")
    n_snips = stack.shape[2]
    track_fs = np.ones(n_snips) * np.nan
    for i in tqdm(range(n_snips)):
        snip = stack[:, :, i]
        fs = np.multiply(kernel, snip)
        track_fs[i] = np.nanmean(fs)
    return track_fs


def generate_scharr_score(stack, kernel=None):
    """
    Generate Scharr score for the stack.
    Convolves the stack with Scharr filter and then convolves it with the filter.
    :param stack: input stack (3D numpy array)
    :param kernel: kernel (or mask, 2D numpy array) for convolution,
                    size should correspond to the first two dimensions of stack
    :return: track with Scharr score (numpy array)
    """
    logger.info(f"Generating Scharr score...")
    n_snips = stack.shape[2]
    track_scharr = np.ones(n_snips) * np.nan
    for i in tqdm(range(n_snips)):
        snip = stack[:, :, i].copy()
        if not kernel is None:
            snip = snip * kernel
        grad = signal.convolve2d(snip, FILTER_SCHARR, boundary="symm", mode="same")
        scharr = np.nanmax(np.absolute(grad))
        track_scharr[i] = scharr
    return track_scharr


def get_peaks_prominence(track):
    """Very simple wrapper to get the prominence of peaks (no filtering).
    Takes 1D numpy array as input."""
    logger.info(f"Finding peaks in track...")
    poss, proms = peaks.find_peak_prominence(track)
    ins_prom_track = np.zeros_like(track) * np.nan
    ins_prom_track[poss] = proms
    return ins_prom_track


#######################
### Fountains analysis:
#######################
def inverted_trinagle_mask(alpha, h, fill_pos=1, fill_neg=0):
    """
    Return square mask with a single fountain (from lower left to upper right).
    Note that obtuse angles cannot be used.
    :param alpha: angle in radians
    :param h:
    :param fill_pos: values for positive regions of the mask
    :param fill_neg: values for negative regions of the mask
    :return: 2D numoy array with a mask. Resulting size will be 2*size + 1.
    """
    mtx = np.zeros([h, h]) + fill_neg
    for x in range(h):
        bgn = np.trunc((x + 1) * np.tan((np.pi / 2 - alpha) / 2)).astype(int)
        end = np.ceil((x + 0.5) / np.tan((np.pi / 2 - alpha) / 2)).astype(int)
        bgn = bgn if bgn > 0 else 0
        end = end if end < h else h
        mtx[x, bgn:end] = fill_pos
    return np.flip(mtx, axis=1).T


def double_triangle_mask(alpha, size, fill_pos=1, fill_neg=0):
    """
    Return square mask with double fountains.
    :param alpha: angle in radians
    :param size: half-size of the square window
    :param fill_pos: values for positive regions of the mask
    :param fill_neg: values for negative regions of the mask
    :return: 2D numoy array with a mask. Resulting size will be 2*size + 1.
    """
    mask_triangle = inverted_trinagle_mask(
        alpha, size, fill_pos=fill_pos, fill_neg=fill_neg
    )
    mask_double_triangle = np.vstack(
        [
            np.hstack([np.ones([size, size]) * (fill_neg), mask_triangle]),
            np.hstack([np.flip(mask_triangle), np.ones([size, size]) * (fill_neg)]),
        ]
    )
    mask = (fill_neg) * np.ones([2 * size + 1, 2 * size + 1])
    mask[np.triu_indices(2 * size + 1, 1)] = mask_double_triangle[
        np.triu_indices(2 * size, 0)
    ]
    mask[np.tril_indices(2 * size + 1, -1)] = mask_double_triangle[
        np.tril_indices(2 * size, 0)
    ]

    return mask


#######################
### Eigenvectors analysis of snips:
#######################
def _get_eig(snip, n_eigs):
    snip = snip.copy()
    snip[np.isnan(snip)] = 0
    eigvecs, eigvals = numutils.get_eig(snip, n_eigs, mask_zero_rows=True)
    eigvecs /= np.sqrt(np.nansum(eigvecs ** 2, axis=1))[:, None]
    eigvecs *= np.sqrt(np.abs(eigvals))[:, None]
    return eigvecs


def _flip_eigs(eigs, vect):
    """ Flipping the sign accoring to the vector sign. """
    return (eigs.T * np.sign(vect)).T


def _get_components(stack, n_eigs):
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


def _flip_sign(eigvects, n_eigs, reference_eigvects):
    n_snips = eigvects.shape[2]
    snip_size = eigvects.shape[1]
    stack_eigvects = np.ones((n_eigs, snip_size, n_snips)) * np.nan
    for i in range(n_snips):
        eigs = eigvects[:, :, i]
        eigs[np.isnan(eigs)] = 0
        if np.sum(eigs) == 0:
            continue
        try:
            vect = np.array(
                [
                    scipy.stats.pearsonr(eigs[i, :], reference_eigvects[i, :])[0]
                    for i in range(n_eigs)
                ]
            )
            stack_eigvects[:, :, i] = flip_eigs(eigs, vect)
        except Exception:
            pass
    return stack_eigvects
