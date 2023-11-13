"""
Rank genes based on various metrics.
"""

import pandas as pd
import logging
from pathlib import Path
from statistics import mean
from typing import Optional
from gene_ranker.dual_dataset import DualDataset
from gene_ranker.ranking_methods import RankingMethod

logging.basicConfig()
log = logging.getLogger("archon")  # An archon was an ancient greek magistrate

def filter_dataset(
    data: pd.DataFrame,
    avg_mean_threshold: Optional[float] = None,
    only_in: Optional[list[str]] = None,
    id_col="gene_id",
) -> pd.DataFrame:
    """Filters a dataset according to simple criteria.

    Filters based on the average expression of the genes (across all samples)
    and a pre-determined list of genes.

    Args:
        data (pd.DataFrame): The dataframe to filter.
        id_col (str): The name of the column with the gene IDs.
        avg_mean_threshold (int or None): Keep only genes with a general average
            expression over this threshold. If None, does not filter.
            Defaults to None.
        only_in (list[str] or None): Keep genes only in this list. If None,
            does not filter. Defaults to None.
    """
    data = data.set_index(id_col)
    if avg_mean_threshold:
        means = data.apply(mean, axis=1)
        keep = means[means > avg_mean_threshold].index
        data = data.filter(items=keep.tolist(), axis=0)

    if only_in:
        data = data.filter(items=only_in, axis=0)

    data[id_col] = data.index
    data = data.reset_index(drop=True)
    
    col_order = [id_col] + [x for x in data.columns if x != id_col]
    data = data[col_order]

    return data


def run_method(
    case_matrix: Path, control_matrix: Path, method: RankingMethod
) -> pd.DataFrame:
    """Run a RankingMethod on two frames.

    Args:
        case_matrix (Path): Path to the case matrix to be read. In `csv` format.
        control_matrix (Path): Same as above, with the control matrix.
        method (RankingMethod): A valid RankingMethod.
    """
    case_matrix: pd.DataFrame = pd.read_csv(case_matrix)
    control_matrix: pd.DataFrame = pd.read_csv(control_matrix)

    dual_dataset = DualDataset(case=case_matrix, control=control_matrix)

    result = method.exec(dual_dataset=dual_dataset)

    return result

