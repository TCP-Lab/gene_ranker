import pytest
from pathlib import Path
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
gene_1,-2.0
gene_2,0.7914728801198858
gene_3,9.129127819482877
"""
expected_bws_test = """\
gene_id,ranking
gene_1,-2.0
gene_2,0.7914728801198858
gene_3,9.129127819482877
"""

# This is before the change of correcting the normalization log2
# gene_1,-2.1213203435596424
# gene_2,0.8142253935739803
# gene_3,1.9986492873450863

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
    args = [case_data_path, control_data_path, "fold_change", "--output-file", target]
    args = [str(x) for x in args]

    bin(args)

    with target.open("r") as stream:
        output = stream.read()

    assert output == expected_fold_change

def test_integration_cohen(tmp_path, case_data_path, control_data_path):
    target = tmp_path / "output.csv"
    args = [case_data_path, control_data_path, "norm_cohen_d", "--output-file", target]
    args = [str(x) for x in args]

    bin(args)

    with target.open("r") as stream:
        output = stream.read()

    assert output == expected_cohen_d

def test_integration_bws(tmp_path, case_data_path, control_data_path):
    target = tmp_path / "output.csv"
    args = [case_data_path, control_data_path, "bws_test", "--output-file", target]
    args = [str(x) for x in args]

    bin(args)

    with target.open("r") as stream:
        output = stream.read()

    assert output == expected_bws_test
