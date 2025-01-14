from io import StringIO
from pathlib import Path

import pandas as pd
import pytest

from gene_ranker.bin import bin

control_data = """\
gene_id,sample_4,sample_5,sample_6
gene_1,6.5,4.0,2.2
gene_2,1.6,0.1,0.1
gene_3,0.0,0.15,0.26
"""
case_data = """\
gene_id,sample_1,sample_2,sample_3
gene_1,2.5,1.0,1.2
gene_2,0.1,0,3.2
gene_3,6.0,3.2,5.01
"""
expected_fold_change = """\
gene_id,ranking
gene_1,-2.666666666666667
gene_2,0.5
gene_3,4.6
"""
expected_cohen_d = """\
gene_id,ranking
gene_1,-1.634014740411334
gene_2,0.35093120317179816
gene_3,4.562438459962434
"""
expected_bws_test = """\
gene_id,ranking
gene_1,-1.44444444444
gene_2,-0.592592592592
gene_3,2.629629629629629
"""
expected_deseq_shrink = """\
gene_id,ranking
gene_1,0.0010467297
gene_2,-0.150855523
gene_3,-8.7566809588
"""


def compare_csvs(one: str, two: str):
    one = pd.read_csv(StringIO(one))
    two = pd.read_csv(StringIO(two))

    pd.testing.assert_frame_equal(one, two, check_like=True, check_exact=False, atol=1e-5)


@pytest.fixture
def case_data_path(tmp_path: Path):
    target = tmp_path / "test_case.csv"
    with target.open("w+") as stream:
        stream.write(case_data)

    return target


@pytest.fixture
def control_data_path(tmp_path: Path):
    target = tmp_path / "test_control.csv"
    with target.open("w+") as stream:
        stream.write(control_data)

    return target


def test_integration_fold_change(tmp_path, case_data_path, control_data_path):
    target = tmp_path / "output.csv"
    args = [case_data_path, control_data_path, "--output-file", target, "fold_change"]
    args = [str(x) for x in args]

    bin(args)

    with target.open("r") as stream:
        output = stream.read()

    compare_csvs(output, expected_fold_change)


def test_integration_cohen(tmp_path, case_data_path, control_data_path):
    target = tmp_path / "output.csv"
    args = [case_data_path, control_data_path, "--output-file", target, "cohen_d"]
    args = [str(x) for x in args]

    bin(args)

    with target.open("r") as stream:
        output = stream.read()

    compare_csvs(output, expected_cohen_d)


def test_integration_bws(tmp_path, case_data_path, control_data_path):
    target = tmp_path / "output.csv"
    args = [case_data_path, control_data_path, "--output-file", target, "bws_test"]
    args = [str(x) for x in args]

    bin(args)

    with target.open("r") as stream:
        output = stream.read()

    compare_csvs(output, expected_bws_test)


def test_integration_deseq(tmp_path, case_data_path, control_data_path):
    target = tmp_path / "output.csv"
    args = [case_data_path, control_data_path, "--output-file", target, "deseq_shrinkage"]
    args = [str(x) for x in args]

    bin(args)

    with target.open("r") as stream:
        output = stream.read()

    compare_csvs(output, expected_deseq_shrink)
