import numpy as np
import pandas as pd
from numpy import ma
from scipy.stats import find_repeats

from gene_ranker.dual_dataset import DualDataset
from gene_ranker.methods.base import fail_if_empty


def rankdata(data):
    def _rank1d(data, use_missing=False):
        n = data.count()
        rk = np.empty(data.size, dtype=float)
        idx = data.argsort()
        rk[idx[:n]] = np.arange(1, n + 1)

        if use_missing:
            rk[idx[n:]] = (n + 1) / 2.0
        else:
            rk[idx[n:]] = 0

        repeats = find_repeats(data.copy())
        for r in repeats[0]:
            condition = (data == r).filled(False)
            rk[condition] = rk[condition].max()
        return rk

    data = ma.array(data, copy=False)
    if data.ndim > 1:
        return _rank1d(data.ravel()).reshape(data.shape)
    else:
        return _rank1d(data)


def bws_score(x, y, alternative):
    axis = 0
    z = rankdata(np.concatenate((x, y)))
    x, y = z[: len(x)], z[len(x) :]

    Ri, Hj = np.sort(x, axis=axis), np.sort(y, axis=axis)
    n, m = Ri.shape[axis], Hj.shape[axis]
    i, j = np.arange(1, n + 1), np.arange(1, m + 1)

    Bx_num = Ri - (m + n) / n * i
    By_num = Hj - (m + n) / m * j

    if alternative == "two-sided":
        Bx_num *= Bx_num
        By_num *= By_num
    else:
        Bx_num *= np.abs(Bx_num)
        By_num *= np.abs(By_num)

    Bx_den = i / (n + 1) * (1 - i / (n + 1)) * m * (m + n) / n
    By_den = j / (m + 1) * (1 - j / (m + 1)) * n * (m + n) / m

    Bx = 1 / n * np.sum(Bx_num / Bx_den, axis=axis)
    By = 1 / m * np.sum(By_num / By_den, axis=axis)

    B = (Bx + By) / 2 if alternative == "two-sided" else (Bx - By) / 2

    return B


@fail_if_empty
def bws_rank(dual_dataset: DualDataset) -> pd.DataFrame:
    dual_dataset.sync()

    case = dual_dataset.case.loc[:, dual_dataset.case.columns != dual_dataset.on]
    control = dual_dataset.control.loc[
        :, dual_dataset.control.columns != dual_dataset.on
    ]
    stats = []

    # The double unpack in because .iterrows() -> (index, row)
    for (_, case), (_, control) in zip(case.iterrows(), control.iterrows()):
        # I needed to get the private method for the statistic that is then
        # resampled many times for the test. The signature therefore is a bit
        # cryptic. The last 0 is the axis.
        stats.append(bws_score(case, control, "one-sided"))

    return pd.DataFrame(
        {dual_dataset.on: dual_dataset.merged[dual_dataset.on], "ranking": stats}
    )
