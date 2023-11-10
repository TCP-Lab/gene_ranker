## Testing
# Just a sprinkle of tests
import pytest
import pandas as pd
from gene_ranker.dual_dataset import DualDataset

@pytest.fixture
def test_case_data():
    return pd.DataFrame(
        {
            "gene_id": ["gene_1", "gene_2", "gene_3"],
            "sample_1": [2.5, 0.1, 6.0],
            "sample_2": [1.0, 0, 3.2],
            "sample_3": [1.2, 0.5, 5.01],
        }
    )


@pytest.fixture
def test_control_data():
    return pd.DataFrame(
        {
            "gene_id": ["gene_1", "gene_2", "gene_3"],
            "sample_4": [6.5, 1.6, 0.0],
            "sample_5": [4.0, 0.1, 0.15],
            "sample_6": [2.2, 0.1, 0.26],
        }
    )

def test_dual_dataset_istantiation(test_case_data, test_control_data):
    dual_data = DualDataset(
        case=test_case_data, control=test_control_data, on="gene_id"
    )
    assert dual_data.case.equals(test_case_data)
    assert dual_data.control.equals(test_control_data)


def test_dual_dataset_merge(test_case_data, test_control_data):
    dual_data = DualDataset(
        case=test_case_data, control=test_control_data, on="gene_id"
    )

    expected_merge = pd.DataFrame(
        {
            "gene_id": ["gene_1", "gene_2", "gene_3"],
            "sample_1": [2.5, 0.1, 6.0],
            "sample_2": [1.0, 0, 3.2],
            "sample_3": [1.2, 0.5, 5.01],
            "sample_4": [6.5, 1.6, 0.0],
            "sample_5": [4.0, 0.1, 0.15],
            "sample_6": [2.2, 0.1, 0.26],
        }
    )

    assert dual_data.merged.equals(expected_merge)


def test_dual_dataset_consistency(test_control_data, test_case_data):
    dual_data = DualDataset(
        case=test_case_data, control=test_control_data, on="gene_id"
    )

    secondary_data = pd.DataFrame(
        {"gene_id": ["gene_1", "gene_2"], "sample_7": [0.0, 0.0]}
    )

    secondary_merge = pd.DataFrame(
        {
            "gene_id": ["gene_1", "gene_2"],
            "sample_7": [0.0, 0.0],
            "sample_4": [6.5, 1.6],
            "sample_5": [4.0, 0.1],
            "sample_6": [2.2, 0.1],
        }
    )

    assert dual_data.case.equals(test_case_data)
    assert dual_data.control.equals(test_control_data)

    dual_data.case = secondary_data

    assert dual_data.case.equals(secondary_data)
    assert dual_data.control.equals(test_control_data)
    assert dual_data.merged.equals(secondary_merge)


def test_dual_dataset_detangle(test_case_data, test_control_data):
    dual_data = DualDataset(
        case=test_case_data, control=test_control_data, on="gene_id"
    )

    expected_merge = pd.DataFrame(
        {
            "gene_id": ["gene_1", "gene_2", "gene_3"],
            "sample_1": [2.5, 0.1, 6.0],
            "sample_2": [1.0, 0, 3.2],
            "sample_3": [1.2, 0.5, 5.01],
            "sample_4": [6.5, 1.6, 0.0],
            "sample_5": [4.0, 0.1, 0.15],
            "sample_6": [2.2, 0.1, 0.26],
        }
    )

    assert dual_data.case.equals(test_case_data)
    assert dual_data.control.equals(test_control_data)
    assert dual_data.merged.equals(expected_merge)
    dual_data._case = None
    dual_data._control = None
    assert dual_data.case.equals(test_case_data)
    assert dual_data.control.equals(test_control_data)

