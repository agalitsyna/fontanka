import click
import click_log
from . import cli

import pandas as pd
import numpy as np

import cooler
from cooltools import snipping
from ..lib import *

# Set up logging:
from . import get_logger
logger = get_logger(__name__)

@cli.command()
@click_log.simple_verbosity_option(logger)
@click.argument("cool_path", metavar="COOL_PATH", type=str)
@click.argument("output_path", metavar="OUTPUT_PATH", type=str)
@click.option(
    "--regions",
    help="File with regions. Required columns: chrom, start, end; tsv format with no header",
    type=str,
)
@click.option(
    "--window-size",
    "-W",
    help="Window size for fountains, in basepairs",
    type=int,
)
@click.option(
    "--angle",
    "-A",
    help="Angle between the sides of the fountain. E.g., np.pi/4 is approximately 0.7854",
    type=np.float32,
)
@click.option(
    "--snips",
    "-S",
    help="npy object with stored whole-genome snips",
    type=str,
    default=None,
)
@click.option(
    "--nthreads",
    "-p",
    help="Number of threads to use for generating snips",
    type=int,
    default=1,
)
@click.option(
    "--store-snips",
    help="File for storing the snips. Will have .npy format.",
    type=str,
    default=None,
)

def call_fountains(
    cool_path, output_path, regions, window_size, angle, snips, nthreads, store_snips
):

    logger.info(
        f"""Running fountain calling for: {cool_path}, 
fountain angle:{angle}, window size: {window_size}"""
    )

    # load cooler
    clr = cooler.Cooler(cool_path)
    bins = clr.bins()[:]
    resolution_bp = clr.binsize

    assert window_size%resolution_bp==0, "Window size should be divisible by resolution"

    # load chromosome regions
    chroms_regions = pd.read_csv(regions)
    chroms_regions.loc[:, "name"] = chroms_regions.apply(
        lambda x: f"{x.chrom}:{x.start}-{x.end}", axis=1
    )

    # create genome-wide windows:
    windows = snipping.make_bin_aligned_windows(
        resolution_bp, bins["chrom"], bins["start"], flank_bp=window_size
    )

    # Assign genomic regions to windows:
    windows = snipping.assign_regions(windows, chroms_regions).reset_index(drop=True)

    # Create and store the stack:
    if not snips is None:
        stack = read_snips(snips)
    else:
        stack = generate_ObsExpSnips(clr, windows, chroms_regions, nthreads=nthreads)
        # Store the stack:
        if not store_snips is None:
            save_snips(stack, store_snips)

    logger.info(f"Finished generating stack, stack shape:{stack.shape} ")

    # Create fountain masks
    # First, calculate the number positive and negative values for normalization:
    n_pos = double_triangle_mask(angle, window_size//resolution_bp, fill_pos=1, fill_neg=0).sum()
    n_neg = double_triangle_mask(angle, window_size//resolution_bp, fill_pos=0, fill_neg=1).sum()
    # All positive values balance out negative values, mask sums up to 0:
    mask_norm = double_triangle_mask(
        angle, window_size//resolution_bp, fill_pos=1 / n_pos, fill_neg=-1 / n_neg
    )
    # Mask sums up to 1:
    mask_pos = double_triangle_mask(angle, window_size//resolution_bp, fill_pos=1 / n_pos, fill_neg=0)

    # Create track with metadata:
    # Fountain score track:
    fs_track = generate_fountain_score(stack, kernel=mask_norm)
    # Scharr score track (for fountain structure only):
    scharr_track = generate_scharr_score(stack, kernel=mask_pos)
    # Scharr score track (for the whole window):
    scharr_track_box = generate_scharr_score(stack)

    # Write the result, note that bad bins are not filtered out:
    metadata = bins[['chrom', 'start', 'end']].copy()
    metadata.loc[:, 'window_start'] = windows['start']
    metadata.loc[:, 'window_end'] = windows['end']
    # Fountain Score (FS) is an average OEE in the fountain divided by average OOE outside of it:
    metadata["FS"] = fs_track
    # Fountain peaks are the prominences of peaks:
    metadata["FS_peaks"] = metadata.groupby("chrom")["FS"].transform(
        get_peaks_prominence
    )
    # Average Scharr in the fountain:
    metadata["Scharr"] = scharr_track
    # Average Scharr gor a box:
    metadata["Scharr_box"] = scharr_track_box

    ######################
    ### Store the dataset
    ######################
    logger.info(f"Writing generated metadata into: {output_path}")
    metadata.to_csv(output_path, sep="\t")
