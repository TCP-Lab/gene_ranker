import pytest
from gene_ranker.ranker import filter_dataset
from gene_ranker.dual_dataset import DualDataset
from gene_ranker.ranking_methods import fold_change_ranking, move_col_to_front
import pandas as pd
from pandas.testing import assert_frame_equal
import pydeseq2
from io import StringIO

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

def test_deseq2_norm():
    # this is a test to see if we get the same results as the R-based deseq2
    # I ran this analysis with this data and got this result.
    # I don't use the test data from the other modules as it has too many
    # zeroes.

    norm_data = """\
"gene_id","sample_1","sample_2","sample_3","sample_4","sample_5","sample_6"
"gene_1",7.70938976055195,2.61532097202366,6.86671430434122,14.6157180877131,7.70938976055195,2.61532097202366
"gene_2",2.56979658685065,2.61532097202366,0.965631699047984,2.56979658685065,2.56979658685065,2.61532097202366
"gene_3",2.56979658685065,2.61532097202366,3.43335715217061,0.803061433390828,2.56979658685065,2.61532097202366
"""
    data = """\
"gene_id","sample_1","sample_2","sample_3","sample_4","sample_5","sample_6"
"gene_1",6,1,64,91,3,1
"gene_2",2,1,9,16,1,1
"gene_3",2,1,32,5,1,1
"""
    metadata = """\
"name","status"
"sample_1","case"
"sample_2","case"
"sample_3","case"
"sample_4","control"
"sample_5","control"
"sample_6","control"
"""
    
    norm_data = pd.read_csv(StringIO(norm_data))
    data = pd.read_csv(StringIO(data)).set_index("gene_id")
    data = data.transpose()

    print(data)

    norm_counts = pydeseq2.preprocessing.deseq2_norm(data)[0].transpose()
    norm_counts["gene_id"] = norm_counts.index
    norm_counts = norm_counts.reset_index(drop=True)
    print(norm_counts)
    print(norm_data)
    assert_frame_equal(norm_counts, norm_data, check_like=True)

