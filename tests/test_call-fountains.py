from click.testing import CliRunner
from fontanka.cli import cli
import os.path as op


def test_call_fountains_cli(request, tmpdir):

    import cooltools

    cool_file = cooltools.download_data(
        "HFF_MicroC", cache=False, data_dir=tmpdir
    )

    regions_file = op.join(tmpdir, "regions.txt")
    with open(regions_file, "w") as f:
        f.write("chrom\tstart\tend\nchr2\t0\t242193529\nchr17\t0\t83257441")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "call-fountains",
            cool_file + "::resolutions/10000",
            "output.fountains.tsv",
            "--regions",
            regions_file,
            "-A",
            0.7854,
            "-p",
            5,
            "-W",
            200_000,
        ],
    )
    assert result.exit_code == 0, result.output
