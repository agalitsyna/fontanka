import click
import click_log
from . import cli

import pandas as pd
import numpy as np
import cooltools.lib.io

import cooler
from cooltools.api import snipping
from ..lib import read_snips, generate_similarity_score, generate_scharr_score, get_peaks_prominence

# Set up logging:
from . import get_logger
logger = get_logger(__name__)

@cli.command()
@click_log.simple_verbosity_option(logger)
@click.argument("cool_path", metavar="COOL_PATH", type=str)
@click.argument("output_path", metavar="OUTPUT_PATH", type=str)
@click.option(
    "--view",
    "--regions",
    help="Path to a BED file which defines which regions of the chromosomes to use. "
    " Note that '--regions' is old-style Open2C name, use '--view' instead. ",
    default=None,
    type=str
)
@click.option(
    "--snips",
    "-S",
    help="NPZ object with stored whole-genome snips",
    type=str,
    required=True,
)
@click.option(
    "--mask",
    "-M",
    help="npy object with stored mask",
    type=str,
    default=None,
)
@click.option(
    "--measure",
    "-m",
    help="Measure for calculating similarity (corr, spearmanr, mse, or mult)",
    type=str,
    default="corr",
)
@click.option(
    "--window-size",
    "-W",
    help="Window size for fountains, in basepairs",
    type=int,
)
def call_mask(
    cool_path, output_path, view, snips, window_size, mask, measure
):

    logger.info(
        f"""Running fountain calling for: {cool_path}, window size: {window_size}
        mask location:{mask}"""
    )

    # load cooler
    clr = cooler.Cooler(cool_path)
    bins = clr.bins()[:]
    resolution_bp = clr.binsize

    assert (
        window_size % resolution_bp == 0
    ), "Window size should be divisible by resolution"

    # load chromosome regions
    chroms_viewframe = cooltools.lib.io.read_viewframe_from_file(view, clr, check_sorting=True)

    # create genome-wide windows:
    windows = snipping.make_bin_aligned_windows(
        resolution_bp, bins["chrom"], bins["start"], flank_bp=window_size
    )

    # Assign genomic view to windows:
    windows = snipping.assign_view_auto(windows, chroms_viewframe).reset_index(drop=True)

    # Store the stack:
    stack = read_snips(snips)

    # Read fountain mask:
    mask = np.load(mask)

    # Create track with metadata:
    # Similarity track:
    FS_track = generate_similarity_score(stack, mask, measure)
    # Scharr score track (for the whole window):
    scharr_track_box = generate_scharr_score(stack)

    # Write the result, note that bad bins are not filtered out:
    metadata = bins[["chrom", "start", "end"]].copy()
    metadata.loc[:, "window_start"] = windows["start"]
    metadata.loc[:, "window_end"] = windows["end"]
    
    # Fountain Score (FS) is an average OEE in the fountain divided by average OOE outside of it:
    metadata["FS"] = FS_track

    # Fountain peaks are the prominences of peaks:
    metadata["FS_peaks"] = metadata.groupby("chrom")["FS"].transform(
        get_peaks_prominence
    )
    # Average Scharr gor a box:
    metadata["Scharr_box"] = scharr_track_box

    ######################
    ### Store the dataset
    ######################
    logger.info(f"Writing generated metadata into: {output_path}")
    metadata.to_csv(output_path, sep="\t")
