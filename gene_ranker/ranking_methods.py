from gene_ranker.dual_dataset import DualDataset
import pandas as pd
from statistics import mean
from dataclasses import dataclass
from typing import Callable, Optional
from shutil import which
from pydeseq2.ds import DeseqDataSet, DeseqStats
from multiprocessing import cpu_count
from math import floor

class MissingExternalDependency(Exception):
    """Raised when an external dependency is missing"""

def check_external_command(command: str, path: None) -> None:
    """Check the presence of a callable for sys-created processes.

    Raises:
        MissingExternalDependency: If the callable is not found.
    """
    if which(command, path=path):
        return
    raise MissingExternalDependency(
        f"External dep '{command}' could not be found in '{path or '$PATH'}'"
    )

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

    print(data)
    print(metadata)

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

    print(shrunk)

    return shrunk


def norm_cohen_d_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    pass


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
