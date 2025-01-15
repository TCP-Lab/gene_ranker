from statistics import mean
import logging

import pandas as pd
import numpy as np
from numpy import std

from gene_ranker.dual_dataset import DualDataset
from gene_ranker.methods.base import fail_if_empty

log = logging.getLogger(__name__)


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

    if np.any(noise == 0):
        log.warn("Some noise values are 0. Setting them to a very small value")
        noise[noise == 0] = 1e-5

    frame = pd.DataFrame(
        {dual_dataset.on: dual_dataset.case[dual_dataset.on], "ranking": signal / noise}
    )

    return frame
