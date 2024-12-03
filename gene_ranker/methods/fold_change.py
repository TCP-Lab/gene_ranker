from statistics import mean

import pandas as pd

from gene_ranker.dual_dataset import DualDataset
from gene_ranker.methods.base import fail_if_empty


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

    frame = pd.DataFrame(
        {dual_dataset.on: dual_dataset.case[dual_dataset.on], "ranking": fcs}
    )

    return frame
