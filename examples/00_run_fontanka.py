import matplotlib as mpl

mpl.use("Agg")
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid import make_axes_locatable

import seaborn as sns
import proplot

from fontanka import *

from cooltools.lib import numutils
import bioframe

import warnings

warnings.filterwarnings("ignore")

labels = [
    "Wild-Type_5.3",
    "NP",
    "PS",
    "SN",
    "sperm",
    "TR",
    "MZnanog_2.75",
    "MZnanog_5.3",
    "MZsox19b_5.3",
    "MZspg_2.75",
    "MZspg_5.3",
    "Wild-Type_11",
    "Wild-Type_25",
    "Wild-Type_2.75",
    "Wild-Type_5.3",
    "WT",
    "embryos_2.3hpf",
    "embryos_24hpf",
    "embryos_4hpf",
    "embryos_8hpf",
]

PATH = "/home/agalicina/DANIO/HIC/data_danrer11/distiller/results_danrer11/coolers_library_group/"

for size, flank_bp, resolution in [
    (20, 200_000, 10_000),
    (25, 250_000, 10_000),
    (50, 500_000, 10_000),
    (10, 250_000, 25_000),
    (10, 250_000, 25_000),
    (40, 1_000_000, 25_000),
    (10, 500_000, 50_000),
    (20, 1_000_000, 50_000),
]:

    for angle_mnemo, angle in [
        ("pi4", np.pi / 4),
        ("pi3", np.pi / 3),
        ("pi6", np.pi / 6),
        ("2pi3", 2 * np.pi / 3),
    ]:

        for label in labels:
            print(label, angle_mnemo, flank_bp, size)

            coolpath = f"{PATH}{label}.danrer11-reduced.mapq_30.1000.mcool::/resolutions/{resolution_bp}"
            clr = cooler.Cooler(coolpath)
            bins = clr.bins()[:]

            chroms_regions = pd.read_csv(
                "/home/agalicina/DANIO/GENOME/danRer11.armsizes.txt"
            )
            chroms_regions.loc[:, "name"] = chroms_regions.apply(
                lambda x: f"{x.chrom}:{x.start}-{x.end}", axis=1
            )

            # Create genome-wide windows:
            windows = snipping.make_bin_aligned_windows(
                resolution_bp, bins["chrom"], bins["start"], flank_bp=flank_bp
            )

            # Assign genomic regions to windows:
            windows = snipping.assign_regions(windows, chroms_regions).reset_index(
                drop=True
            )

            # Create and store the stack:
            try:
                stack = read_snips(
                    f"data/Snips_WholeGenome_{label}_Res{resolution_bp}_Flank{flank_bp}_Size{size}_Angle-{angle_mnemo}.npy"
                )
            except Exception as e:
                stack = generate_ObsExpSnips(clr, windows, chroms_regions, nthreads=10)
                # Store the stack:
                save_snips(
                    stack,
                    f"data/Snips_WholeGenome_{label}_Res{resolution_bp}_Flank{flank_bp}_Size{size}_Angle-{angle_mnemo}.npy",
                )

            print("Stack size: ", stack.shape)

            # Create fountain masks
            # First, calculate the number positive and negative values for normalization:
            n_pos = double_triangle_mask(angle, size, fill_pos=1, fill_neg=0).sum()
            n_neg = double_triangle_mask(angle, size, fill_pos=0, fill_neg=1).sum()
            # All positive values balance out negative values, mask sums up to 0:
            mask_norm = double_triangle_mask(
                angle, size, fill_pos=1 / n_pos, fill_neg=-1 / n_neg
            )
            # Mask sums up to 1:
            mask_pos = double_triangle_mask(angle, size, fill_pos=1 / n_pos, fill_neg=0)

            # Create track with metadata:
            # Fountain score track:
            fs_track = generate_fountain_score(stack, kernel=mask_norm)
            # Scharr score track (for fountain structure only):
            scharr_track = generate_scharr_score(stack, kernel=mask_pos)
            # Scharr score track (for the whole window):
            scharr_track_box = generate_scharr_score(stack)

            # Write the result, note that bad bins are not filtered out:
            metadata = bins[["chrom", "start", "end"]].copy()
            metadata.loc[:, "window_start"] = windows["start"]
            metadata.loc[:, "window_end"] = windows["end"]
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
            metadata.to_csv(
                f"data/Table_{label}_Res{resolution_bp}_Flank{flank_bp}_Size{size}_Angle-{angle_mnemo}.tsv",
                sep="\t",
            )

            #########################
            ### Very simple filtering
            #########################
            FS_threshold = np.nanpercentile(metadata["FS_peaks"], 80)

            #################
            ### Visualization
            #################
            # Figure 1:
            snip = stack[:, :, 460]

            f, axs = plt.subplots(
                figsize=(10, 6), nrows=1, ncols=2, sharex=True, sharey=True
            )

            ax = axs[0]
            ax.set_title("Snip")
            im = ax.matshow(np.log2(snip), cmap="RdBu_r", vmax=1.5, vmin=-1.5)

            ax = axs[1]
            ax.set_title("Fountain mask")
            im1 = ax.matshow(mask_norm, cmap="fall")

            ticks = np.arange(0, len(mask_norm), 10)
            ticklabels = (ticks * resolution_bp - flank_bp) // 1000
            for ax in axs:
                ax.set_yticks(ticks)
                ax.set_yticklabels(ticklabels)
                ax.set_xticks(ticks)
                ax.set_xticklabels([])
                ax.xaxis.tick_bottom()

            f.tight_layout()
            f.savefig(
                f"data/Snip_408_{label}_Res{resolution_bp}_Flank{flank_bp}_Size{size}_Angle-{angle_mnemo}.png"
            )

            # Figure 2:
            mtx = clr.matrix(as_pixels=False, balance=True).fetch("chr1")
            for i in range(2):
                np.fill_diagonal(mtx[i:, :], np.nan)
                np.fill_diagonal(mtx[:, i:], np.nan)
            image = np.log2(
                numutils.observed_over_expected(mtx, mask=np.isfinite(mtx))[0]
            )

            # Set up the parameters for plotting:
            start, end = 0, 500
            xticks = np.arange(start, end + 1, 100)
            xticks_kb = (xticks * resolution_bp) / 1000000

            f, ax = plt.subplots(
                figsize=(10, 15),
                sharex=True,
            )

            im = ax.matshow(
                image[start:end, start:end], vmax=1.5, vmin=-1.5, cmap="RdBu_r"
            )
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            plt.colorbar(im, cax=cax, label="corrected frequencies")
            ax.set_title("Observed Over Expected Matrix")
            ax.xaxis.set_visible(False)

            ax1 = divider.append_axes("bottom", size="25%", pad=0.1, sharex=ax)
            ax1.plot(metadata["cis_coverage"], label="cis")
            ax1.plot(metadata["total_coverage"], label="total")
            ax1.set_xlim([start, end])
            ax1.set_ylim([0, np.nanpercentile(metadata["total_coverage"], 99)])
            ax1.set_ylabel("Coverage")
            ax1.legend()

            ax2 = divider.append_axes("bottom", size="10%", pad=0.1, sharex=ax)
            ax2.plot(metadata["cis_coverage"] / metadata["total_coverage"])
            ax2.set_xlim([start, end])
            ax2.set_ylabel("coverage ratio")

            ax3 = divider.append_axes("bottom", size="25%", pad=0.1, sharex=ax)
            ax3.plot(metadata["FS"])
            ax3.set_xlim([start, end])
            ax3.set_ylim(
                [
                    np.nanpercentile(metadata["FS"], 1),
                    np.nanpercentile(metadata["FS"], 99),
                ]
            )
            ax3.set_ylabel("Fountain score")

            for i in range(start, end):
                fs = metadata["FS_peaks"][i]
                if not np.isnan(fs) and fs > FS_threshold:
                    ax3.plot(
                        i - start, metadata["FS"][i], marker="v", color="blue", ms=3
                    )

            ax4 = divider.append_axes("bottom", size="25%", pad=0.1, sharex=ax)
            ax4.plot(metadata["Scharr"])
            ax4.set_xlim([start, end])
            ax4.set_ylim(
                [
                    np.nanpercentile(metadata["Scharr"], 1),
                    np.nanpercentile(metadata["Scharr"], 99),
                ]
            )
            ax4.set_ylabel("Scharr score")
            ax4.set(xticks=xticks, xticklabels=xticks_kb, xlabel="Mb")

            f.savefig(
                f"data/Region_chr1_{label}_Res{resolution_bp}_Flank{flank_bp}_Size{size}_Angle-{angle_mnemo}.png"
            )

            # Fountain score selection and avfountain:
            selected = metadata.dropna(axis=0, subset=["FS_peaks"]).query(
                f"FS_peaks>{FS_threshold}"
            )

            avg_fountain = np.nanmean(stack[:, :, selected.index], axis=2)
            eigs_fountain = get_eig(avg_fountain, n_eigs=snip_size)
            # Phasing the eigenvectors by the first value, which it arbitrary
            vect = eigs_fountain[:, 0]
            eigs_fountain = flip_eigs(eigs_fountain, vect)

            # Figure 3:
            f, axs = plt.subplots(
                figsize=(10, 6), nrows=1, ncols=2, sharex=True, sharey=True
            )

            ax = axs[0]
            ax.set_title(f"Average fountain, N={len(selected)}")
            im = ax.matshow(np.log2(avg_fountain), cmap="RdBu_r", vmax=1.5, vmin=-1.5)

            ax = axs[1]
            ax.set_title("eigenvectors")
            im1 = ax.matshow(eigs_fountain, cmap="RdBu")

            ticks = np.arange(0, snip_size, 10)
            ticklabels = (ticks * resolution_bp - flank_bp) // 1000
            for ax in axs:
                ax.set_yticks(ticks)
                ax.set_yticklabels(ticklabels)
                ax.set_xticks(ticks)
                ax.set_xticklabels([])
                ax.xaxis.tick_bottom()

            f.tight_layout()
            f.savefig(
                f"data/AvFountain_{label}_Res{resolution_bp}_Flank{flank_bp}_Size{size}_Angle-{angle_mnemo}.png"
            )

            # Figure 4:
            for n in [460]:
                snip = stack[:, :, n]
                eigs = get_eig(snip, n_eigs=snip_size)
                # Phasing the eigenvectors by sign of pearson correlation with idealized fountain:
                vect = np.array(
                    [
                        scipy.stats.pearsonr(eigs_fountain[i, :], eigs[i, :])[0]
                        for i in range(snip_size)
                    ]
                )
                eigs = flip_eigs(eigs, vect)

                f, axs = plt.subplots(
                    figsize=(10, 6), nrows=1, ncols=2, sharex=True, sharey=True
                )

                ax = axs[0]
                ax.set_title("Random fountain")
                im = ax.matshow(np.log2(snip), cmap="RdBu_r", vmax=1.5, vmin=-1.5)

                ax = axs[1]
                ax.set_title("eigenvectors")
                im1 = ax.matshow(eigs, cmap="RdBu")

                ticks = np.arange(0, snip_size, 10)
                ticklabels = (ticks * resolution_bp - flank_bp) // 1000
                for ax in axs:
                    ax.set_yticks(ticks)
                    ax.set_yticklabels(ticklabels)
                    ax.set_xticks(ticks)
                    ax.set_xticklabels([])
                    ax.xaxis.tick_bottom()

                f.tight_layout()
                f.savefig(
                    f"data/Snip_{n}-witheig_{label}_Res{resolution_bp}_Flank{flank_bp}_Size{size}_Angle-{angle_mnemo}.png"
                )
