import os
import os.path as op
import numpy as np
import pytest
from click.testing import CliRunner

from fontanka.cli import cli
from fontanka.lib.utils import double_triangle_mask

import cooltools
import subprocess


@pytest.fixture(scope="session")
def hff_cool(tmp_path_factory):
    """Download a small public cooler for exercising CLI commands."""
    data_dir = tmp_path_factory.mktemp("data")
    return cooltools.download_data(
        "HFF_MicroC", cache=False, data_dir=data_dir
    )


@pytest.fixture(scope="session")
def regions_file(tmp_path_factory):
    """Minimal viewframe covering two chromosomes."""
    path = tmp_path_factory.mktemp("regions") / "regions.txt"
    path.write_text("chr2\t0\t1000000\nchr17\t1000000\t3000000\n")
    return path


@pytest.fixture(scope="session")
def window_size():
    return 200_000


@pytest.fixture(scope="session")
def stack_file(tmp_path_factory, hff_cool, regions_file, window_size):
    """Generate a stack via slice-windows to reuse in other CLI tests."""
    runner = CliRunner()
    out = tmp_path_factory.mktemp("stacks") / "stack.npy"
    result = runner.invoke(
        cli,
        [
            "slice-windows",
            f"{hff_cool}::resolutions/10000",
            str(out),
            "--regions",
            str(regions_file),
            "-W",
            window_size,
            "-p",
            1,
        ],
    )
    assert result.exit_code == 0, result.output
    assert out.exists()
    return out


@pytest.fixture(scope="session")
def mask_file(tmp_path_factory, stack_file):
    """Build a simple fountain mask sized to the generated snips."""
    stack = np.load(stack_file)
    # stack shape: (n, win, win); derive half-size from window dimension
    half_size = (stack.shape[1] - 1) // 2
    mask = double_triangle_mask(np.pi / 4, half_size, fill_pos=1, fill_neg=0)
    path = tmp_path_factory.mktemp("mask") / "mask.npy"
    np.save(path, mask)
    return path


@pytest.fixture(scope="session")
def expected_file(tmp_path_factory, hff_cool, regions_file):
    """
    Precompute expected using cooltools expected-cis CLI
    (balanced, ignore-diags 2, no-smoothing).
    """
    out = tmp_path_factory.mktemp("expected") / "hff.expected.tsv"

    cmd = [
        "cooltools",
        "expected-cis",
        f"{hff_cool}::resolutions/10000",
        "--regions",
        str(regions_file),
        "-p",
        "1",
        "--clr-weight-name",
        "weight",
        "--ignore-diags",
        "2",
        "-o",
        str(out),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.skip(f"cooltools expected-cis failed: {result.stderr}")
    return out


@pytest.fixture(scope="session")
def expected_file_wholegenome(tmp_path_factory, hff_cool):
    """
    Precompute expected using cooltools expected-cis CLI
    (balanced, ignore-diags 2, no-smoothing).
    """
    out = tmp_path_factory.mktemp("expected") / "hff.expected.wholegenome.tsv"

    cmd = [
        "cooltools",
        "expected-cis",
        f"{hff_cool}::resolutions/10000",
        "-p",
        "1",
        "--clr-weight-name",
        "weight",
        "--ignore-diags",
        "2",
        "-o",
        str(out),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.skip(f"cooltools expected-cis failed: {result.stderr}")
    return out


def test_slice_windows_cli(stack_file):
    """slice-windows should run and emit a stack file."""
    assert stack_file.exists()
    stack = np.load(stack_file)
    assert stack.ndim == 3
    assert stack.shape[1] == stack.shape[2]


def test_apply_fountain_mask_cli(
    hff_cool, stack_file, regions_file, mask_file, window_size
):
    runner = CliRunner()
    out = op.join(op.dirname(stack_file), "fountain_mask.tsv")
    result = runner.invoke(
        cli,
        [
            "apply-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "--regions",
            str(regions_file),
            "-W",
            window_size,
            "--mask",
            str(mask_file),
            "--snips",
            str(stack_file),
            "-p",
            1,
            "--measure",
            "corr",
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_apply_fountain_mask_cli_noregions(
    hff_cool, stack_file, mask_file, window_size
):
    runner = CliRunner()
    out = op.join(op.dirname(stack_file), "fountain_mask_noregions.tsv")
    result = runner.invoke(
        cli,
        [
            "apply-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "-W",
            window_size,
            "--mask",
            str(mask_file),
            "--snips",
            str(stack_file),
            "-p",
            1,
            "--measure",
            "corr",
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_apply_fountain_mask_cli_nostack(
    hff_cool, stack_file, mask_file, window_size
):
    runner = CliRunner()
    out = op.join(op.dirname(stack_file), "fountain_mask_nostack.tsv")
    result = runner.invoke(
        cli,
        [
            "apply-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "-W",
            window_size,
            "--mask",
            str(mask_file),
            "-p",
            1,
            "--measure",
            "corr",
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_apply_fountain_mask_cli_nostack_with_expected(
    hff_cool, stack_file, mask_file, regions_file, expected_file, window_size
):
    runner = CliRunner()
    out = op.join(
        op.dirname(stack_file), "fountain_mask_nostack_with_expected.tsv"
    )
    result = runner.invoke(
        cli,
        [
            "apply-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "-W",
            window_size,
            "--mask",
            str(mask_file),
            "--regions",
            str(regions_file),
            "-p",
            1,
            "--expected",
            str(expected_file),
            "--measure",
            "corr",
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_apply_fountain_mask_cli_nostack_with_expected_wholegenome(
    hff_cool, stack_file, mask_file, expected_file_wholegenome, window_size
):
    runner = CliRunner()
    out = op.join(
        op.dirname(stack_file), "fountain_mask_with_expected_wholegenome.tsv"
    )
    result = runner.invoke(
        cli,
        [
            "apply-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "-W",
            window_size,
            "--mask",
            str(mask_file),
            "-p",
            1,
            "--expected",
            str(expected_file_wholegenome),
            "--measure",
            "corr",
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_apply_binary_fountain_mask_cli(
    hff_cool, stack_file, regions_file, window_size
):
    runner = CliRunner()
    out = op.join(op.dirname(stack_file), "binary_mask.tsv")
    result = runner.invoke(
        cli,
        [
            "apply-binary-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "--regions",
            str(regions_file),
            "-W",
            window_size,
            "-A",
            0.7854,
            "--snips",
            str(stack_file),
            "-p",
            1,
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_apply_binary_fountain_mask_cli_noregions(
    hff_cool, stack_file, window_size
):
    runner = CliRunner()
    out = op.join(op.dirname(stack_file), "binary_mask_noregions.tsv")
    result = runner.invoke(
        cli,
        [
            "apply-binary-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "-W",
            window_size,
            "-A",
            0.7854,
            "--snips",
            str(stack_file),
            "-p",
            1,
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_apply_binary_fountain_mask_cli_nostack(
    hff_cool, stack_file, window_size
):
    runner = CliRunner()
    out = op.join(op.dirname(stack_file), "binary_mask_nostack.tsv")
    result = runner.invoke(
        cli,
        [
            "apply-binary-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "-W",
            window_size,
            "-A",
            0.7854,
            "-p",
            1,
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_apply_binary_fountain_mask_cli_nostack_with_expected(
    hff_cool, stack_file, regions_file, expected_file, window_size
):
    runner = CliRunner()
    out = op.join(
        op.dirname(stack_file), "binary_mask_nostack_with_expected.tsv"
    )
    result = runner.invoke(
        cli,
        [
            "apply-binary-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "--regions",
            str(regions_file),
            "-W",
            window_size,
            "-A",
            0.7854,
            "--expected",
            str(expected_file),
            "-p",
            1,
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_apply_binary_fountain_mask_cli_nostack_with_expected_wholegenome(
    hff_cool, stack_file, expected_file_wholegenome, window_size
):
    runner = CliRunner()
    out = op.join(
        op.dirname(stack_file),
        "binary_mask_nostack_with_expected_wholegenome.tsv",
    )
    result = runner.invoke(
        cli,
        [
            "apply-binary-fountain-mask",
            f"{hff_cool}::resolutions/10000",
            out,
            "-W",
            window_size,
            "-A",
            0.7854,
            "--expected",
            str(expected_file_wholegenome),
            "-p",
            1,
        ],
    )
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)
