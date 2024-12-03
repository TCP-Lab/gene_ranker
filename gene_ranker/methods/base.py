from argparse import ArgumentParser
from dataclasses import dataclass
from functools import wraps
from typing import Callable, Optional

import pandas as pd
from numpy import log2
from pydeseq2.preprocessing import deseq2_norm


def move_col_to_front(data: pd.DataFrame, col_name) -> pd.DataFrame:
    """Move a column to the first position, useful when saving data to disk."""
    cols = data.columns.tolist()
    assert col_name in cols
    cols.insert(0, cols.pop(cols.index(col_name)))
    data = data.loc[:, cols]

    return data


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

    def __post_init__(self):
        if self.parser is None:
            # Set a dummy parser with no options.
            self.parser = ArgumentParser(self.name, description=self.desc)


def norm_with_deseq(data: pd.DataFrame, id_col=None):
    """Normalize a dataframe with the "mean of ratios" method as used by Deseq

    Args:
        data (pandas.DataFrame): A dataframe to normalize.
        id_col (Optional[str]): Optionally, the column with IDs. If not passed,
            assumes the ID column is not present.

    Returns:
        A pandas.DataFrame with normalized counts. The ID column is untouched.
    """
    # Move the ID col to the index if needed
    assert id_col in data.columns
    if id_col:
        data = data.set_index(id_col)
    # Convert back to counts
    data = round((2**data) - 1)
    data = data.map(int)
    data = data.transpose()
    data = deseq2_norm(data)[0].transpose()

    # Return to logged values
    data = log2(data + 1)

    if id_col:
        data[id_col] = data.index
        data.reset_index(drop=True)

    return data


def norm_wrapper(exec: Callable):
    """
    Wrap a standard ranking method to call it with normalized data

    Uses the 'norm_with_deseq' function to normalize the data and then call
    the wrapped method.
    """

    def wrapped(dual_dataset, *args, **kwargs):
        dual_dataset.sync()
        dual_dataset.merged = norm_with_deseq(dual_dataset.merged, dual_dataset.on)

        return exec(dual_dataset, *args, **kwargs)

    return wrapped
