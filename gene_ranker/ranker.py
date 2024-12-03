"""
Rank genes based on various metrics.
"""

import logging
from pathlib import Path
from statistics import mean
from typing import Optional

import pandas as pd

from gene_ranker.dual_dataset import DualDataset
from gene_ranker.methods.base import RankingMethod

log = logging.getLogger(__name__)


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
    case_matrix: Path,
    control_matrix: Path,
    method: RankingMethod,
    shared_col: str = "gene_id",
    extra_args: Optional[dict] = None,
) -> pd.DataFrame:
    """Run a RankingMethod on two frames.

    Args:
        case_matrix (Path): Path to the case matrix to be read. In `csv` format.
        control_matrix (Path): Same as above, with the control matrix.
        method (RankingMethod): A valid RankingMethod.
    """
    case_matrix_data: pd.DataFrame = pd.read_csv(case_matrix)
    control_matrix_data: pd.DataFrame = pd.read_csv(control_matrix)

    log.info(
        f"Loaded a {case_matrix_data.shape[1]} col by {case_matrix_data.shape[0]} rows case matrix from {case_matrix}"
    )
    log.info(
        f"Loaded a {control_matrix_data.shape[1]} col by {control_matrix_data.shape[0]} rows control matrix from {control_matrix}"
    )

    dual_dataset = DualDataset(
        case=case_matrix_data, control=control_matrix_data, on=shared_col
    )

    result = method.exec(dual_dataset=dual_dataset, **extra_args)

    return result
