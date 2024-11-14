from gene_ranker.dual_dataset import DualDataset
import pandas as pd
from statistics import mean
from dataclasses import dataclass
from typing import Callable, Optional
import shutil
from pydeseq2.ds import DeseqDataSet, DeseqStats
import tempfile
import subprocess
from pydeseq2.preprocessing import deseq2_norm
from numpy import log2, std
from scipy.stats._bws_test import _bws_statistic
from functools import wraps


def fail_if_empty(func):
    """Fail fast if case or control are empty

    A decorator for functions that take `dual_dataset`s as input and should
    fail if either case or control matrices are empty.
    """
    @wraps(func)
    def wrap(dual_dataset, *args, **kwargs):    
        case = dual_dataset.case.set_index(dual_dataset.on)
        control = dual_dataset.control.set_index(dual_dataset.on)
        if case.empty:
            raise ValueError("Case matrix is empty. Cannot compute fold change.")
        if control.empty:
            raise ValueError("Control matrix is empty. Cannot compute fold change.")

        return func(dual_dataset, *args, **kwargs)

    return wrap


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

@fail_if_empty
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

    case = dual_dataset.case.set_index(dual_dataset.on)
    control = dual_dataset.control.set_index(dual_dataset.on)
    case_means = case.apply(mean, axis=1)
    control_means = control.apply(mean, axis=1)

    # Assume that the values are logged
    fcs = case_means.to_numpy() - control_means.to_numpy()

    frame = pd.DataFrame({dual_dataset.on: dual_dataset.case[dual_dataset.on], "ranking": fcs})

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
    
    # Return to logged values
    data = log2(data + 1)

    if id_col:
        data[id_col] = data.index
        data.reset_index(drop=True)

    return data

def norm_wrapper(exec: Callable, *args, **kwargs):
    def wrapped(dual_dataset):
        dual_dataset.sync()
        dual_dataset.merged = norm_with_deseq(dual_dataset.merged, dual_dataset.on)

        return exec(dual_dataset, *args, **kwargs)
    return wrapped

@fail_if_empty
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
        quiet=True
    )
    
    data.deseq2()
    stats = DeseqStats(data)
    stats.summary()
    stats.lfc_shrink(coeff="status_case_vs_control")

    shrunk = stats.LFC

    result = pd.DataFrame({
        dual_dataset.on: shrunk.index,
        "ranking": shrunk["status_case_vs_control"]
    })

    return result

@fail_if_empty
def cohen_d_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    if not shutil.which("fast-cohen"):
        raise MissingExternalDependency((
            "Missing the fast-cohen executable."
            " See https://github.com/mrhedmad/fast-cohen/ to download."
        ))

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

@fail_if_empty
def signal_to_noise_ratio(dual_dataset: DualDataset) -> pd.DataFrame:
    dual_dataset.sync()

    case = dual_dataset.case.set_index(dual_dataset.on)
    control = dual_dataset.control.set_index(dual_dataset.on)

    case_means = case.apply(mean, axis=1)
    case_stdev = case.apply(lambda x: std(x, ddof=1), axis=1)
    control_means = control.apply(mean, axis=1)
    control_stdev = control.apply(lambda x: std(x, ddof=1), axis=1)

    # Assume that the values are logged
    signal = case_means.to_numpy() - control_means.to_numpy()
    noise = case_stdev.to_numpy() + control_stdev.to_numpy()

    frame = pd.DataFrame(
        {
            dual_dataset.on: dual_dataset.case[dual_dataset.on],
            "ranking": signal / noise
        }
    )
    
    return frame

@fail_if_empty
def bws_rank(dual_dataset: DualDataset) -> pd.DataFrame:
    dual_dataset.sync()
    
    case = dual_dataset.case.loc[:, dual_dataset.case.columns != dual_dataset.on]
    control = dual_dataset.control.loc[:, dual_dataset.control.columns != dual_dataset.on]
    stats = []

    # The double unpack in because .iterrows() -> (index, row)
    for (_, case), (_, control) in zip(case.iterrows(), control.iterrows()):
        # I needed to get the private method for the statistic that is then
        # resampled many times for the test. The signature therefore is a bit
        # cryptic. The last 0 is the axis.
        stats.append(_bws_statistic(case, control, 'one-sided', 0))

    return pd.DataFrame({
        dual_dataset.on: dual_dataset.merged[dual_dataset.on],
        "ranking": stats
    })


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
        desc = "Use DESeq2-shrunk fold changes. Always normalizes the input"
    ),
    "cohen_d": RankingMethod(
        name = "Cohen's d",
        exec = cohen_d_ranking,
        parser = None,
        desc = "Use the Cohen's d metric"
    ),
    "norm_cohen_d": RankingMethod(
        name = "Normalized Cohen's d",
        exec = norm_wrapper(cohen_d_ranking),
        parser = None,
        desc = "Use a DESeq2-normalized Cohen's d metric"
    ),
    "norm_fold_change": RankingMethod(
        name = "Normalized Fold Change",
        exec = norm_wrapper(fold_change_ranking),
        parser = None,
        desc = "Use a DESeq2-normalized fold change metric"
    ),
    "s2n_ratio": RankingMethod(
        name = "Signal to noise ratio",
        exec = signal_to_noise_ratio,
        parser=None,
        desc = "Use the signal to noise ratio (diff of means divided by variance)"
    ),
    "norm_s2n_ratio": RankingMethod(
        name = "Normalized signal to noise ratio",
        exec = norm_wrapper(signal_to_noise_ratio),
        parser=None,
        desc = "Use the signal to noise ratio metric on normalized data"
    ),
    "bws_test": RankingMethod(
        name = "Baumgartner-Weiss-Schindler test statistic",
        exec = bws_rank,
        parser=None,
        desc = "Use the BWS test statistic, which works well with high N samples"
    ),
    "norm_bws_test": RankingMethod(
        name = "Normalized Baumgartner-Weiss-Schindler test statistic",
        exec = norm_wrapper(bws_rank),
        parser=None,
        desc = "Same as BWS, but on normalized data"
    ),
}

