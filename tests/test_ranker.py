import pytest
from gene_ranker.ranker import filter_dataset
from gene_ranker.dual_dataset import DualDataset
from gene_ranker.ranking_methods import fold_change_ranking
import pandas as pd
from pandas.testing import assert_frame_equal

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

def test_data_filtering_avg(test_case_data):
    new_data: pd.DataFrame = filter_dataset(test_case_data, avg_mean_threshold=1)

    expected = pd.DataFrame(
        {
            "gene_id": ["gene_1", "gene_3"],
            "sample_1": [2.5, 6.0],
            "sample_2": [1.0, 3.2],
            "sample_3": [1.2, 5.01],
        }
    )
    
    assert new_data.equals(expected)

def test_data_filtering_only_in(test_case_data):
    new_data: pd.DataFrame = filter_dataset(test_case_data, only_in=["gene_2", "gene_3"])

    expected = pd.DataFrame(
        {
            "gene_id": ["gene_2", "gene_3"],
            "sample_1": [0.1, 6.0],
            "sample_2": [0, 3.2],
            "sample_3": [0.5, 5.01],
        }
    )

    print(new_data)
    print(expected)
    
    assert new_data.equals(expected)


def test_fold_change_ranking(test_case_data, test_control_data):
    dual_data = DualDataset(case = test_case_data, control=test_control_data)

    result = fold_change_ranking(dual_data)

    expected = pd.DataFrame({
        "gene_id": ["gene_1", "gene_2", "gene_3"],
        "ranking": [-2.66666666, -0.4, 4.6],
    })

    print(result)
    print(expected)
    
    # The `assert` is in the function itself. Used for tolerance of estimates
    # by default allows a tolerance of 1e-5
    assert_frame_equal(result, expected)

