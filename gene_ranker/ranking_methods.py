from gene_ranker.dual_dataset import DualDataset
import pandas as pd
from statistics import mean
from dataclasses import dataclass
from typing import Callable, Optional
import shutil
from pydeseq2.ds import DeseqDataSet, DeseqStats
from multiprocessing import cpu_count
import tempfile
import subprocess
from pydeseq2.preprocessing import deseq2_norm

class MissingExternalDependency(Exception):
    """Raised when an external dependency is missing"""

@dataclass
class RankingMethod:
    """Represents a standard RankingMethod"""
    name: str
    """The human-friendly name of the method"""
    exec: Callable
    """The callable to call with this method"""
    parser: Optional[Callable]
    """An ArgumentParser to use to add options to the callable for this method."""
    desc: Optional[str] = None
    """A human-friendly description of the method."""

def fold_change_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    """Perform a simple fold-change ranking metric.

    Expects log(x + 1) count data, or generally logged data as input.
    For each gene, we compute `mean(case) - mean(control)`, and the result is
    saved.

    Args:
        dual_dataset (DualDataset): A DualDataset to calculate the result from.
    Returns:
        A pd.DataFrame with two columns, a `gene_id` column and a `ranking` column.
    """
    dual_dataset.sync()

    case = dual_dataset.case.set_index("gene_id")
    control = dual_dataset.control.set_index("gene_id")
    case_means = case.apply(mean, axis=1)
    control_means = control.apply(mean, axis=1)

    # Assume that the values are logged
    fcs = case_means.to_numpy() - control_means.to_numpy()

    frame = pd.DataFrame({"gene_id": dual_dataset.case["gene_id"], "ranking": fcs})

    return frame

def move_col_to_front(data: pd.DataFrame, col_name) -> pd.DataFrame:
    cols = data.columns.tolist()
    assert col_name in cols
    cols.insert(0, cols.pop(cols.index(col_name)))
    data = data.loc[:, cols]

    return data

def norm_with_deseq(data: pd.DataFrame, id_col = None):
    # Move the ID col to the index if needed
    assert id_col in data.columns
    if id_col:
        data = data.set_index(id_col)
    # Convert back to counts
    data = round((2 ** data) - 1)
    data = data.map(int)
    data = data.transpose()
    data = deseq2_norm(data)[0].transpose()

    if id_col:
        data[id_col] = data.index
        data.reset_index(drop=True)

    return data

def deseq_shrinkage_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    case_cols = [x for x in dual_dataset.case.columns if x != dual_dataset.on]
    ctrl_cols = [x for x in dual_dataset.control.columns if x != dual_dataset.on]
    labels = ["case"] * len(case_cols) + ["control"] * len(ctrl_cols)
    metadata = pd.DataFrame({
        "sample": case_cols + ctrl_cols,
        "status": labels
    })
    metadata = metadata.set_index("sample")
    data = dual_dataset.merged.set_index(dual_dataset.on)

    data = round((2 ** data) - 1)
    data = data.applymap(int)
    data = data.transpose()

    data = DeseqDataSet(
        counts=data,
        metadata=metadata,
        design_factors="status",
        ref_level=["status", "control"],
        n_cpus=cpu_count(),
        quiet=True
    )
    
    data.deseq2()
    stats = DeseqStats(data)
    stats.run_wald_test()
    stats.lfc_shrink(coeff="status_case_vs_control")

    shrunk = stats.LFC

    result = pd.DataFrame({
        dual_dataset.on: shrunk.index,
        "ranking": shrunk["status_case_vs_control"]
    })

    return result


def norm_cohen_d_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    if not shutil.which("fast-cohen"):
        raise MissingExternalDependency((
            "Missing the fast-cohen executable."
            " See https://github.com/mrhedmad/fast-cohen/ to download."
        ))

    dual_dataset.merged = norm_with_deseq(dual_dataset.merged, dual_dataset.on)

    with tempfile.NamedTemporaryFile() as case, \
            tempfile.NamedTemporaryFile() as control, \
            tempfile.NamedTemporaryFile() as result:
        move_col_to_front(dual_dataset.case, dual_dataset.on).to_csv(case, index=False)
        move_col_to_front(dual_dataset.control, dual_dataset.on).to_csv(control, index=False)

        subprocess.run(["fast-cohen", case.name, control.name, result.name]) 
    
        data = pd.read_csv(result)
        # Rename the 'row_names' col
        data = data.rename(columns = {"row_names": dual_dataset.on, "cohen_d": "ranking"})

    return data


def norm_hedges_g_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    pass


# Tried with an enum but it's just too much of a bother to implement.
# A simple dict is fine, ultimately.
RANKING_METHODS = {
    "fold_change": RankingMethod(
        name = "Fold Change",
        exec = fold_change_ranking,
        parser = None,
        desc = "Use a non-normalized, raw fold change metric."
    ),
    "deseq_shrinkage": RankingMethod(
        name = "DESeq2 Shrinkage",
        exec = deseq_shrinkage_ranking,
        parser = None,
        desc = "Use DESeq2-shrunk fold changes."
    ),
    "norm_cohen_d": RankingMethod(
        name = "Normalized Cohen's D",
        exec = norm_cohen_d_ranking,
        parser = None,
        desc = "Use a DESeq2-normalized Cohen's D metric"
    ),
    "norm_hedges_g": RankingMethod(
        name = "Normalized Hedges' G",
        exec = norm_hedges_g_ranking,
        parser = None,
        desc = "Use a DESeq2-normalized Hedges' G metric"
    ),
}
