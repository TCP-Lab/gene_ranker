from gene_ranker.dual_dataset import DualDataset
import pandas as pd
from statistics import mean
from dataclasses import dataclass
from typing import Callable, Optional
from shutil import which

class MissingExternalDependency(Exception):
    """Raised when an external dependency is missing"""

def check_external_command(command: str, path: None) -> None:
    if which(command, path=path):
        return
    raise MissingExternalDependency(
        f"External dep '{command}' could not be found in '{path or '$PATH'}'"
    )

@dataclass
class RankingMethod:
    name: str
    exec: Callable
    parser: Optional[Callable]
    desc: Optional[str] = None

def fold_change_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
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
    pass


def norm_cohen_d_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    pass


def norm_hedges_g_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    pass


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

