import pandas as pd

from gene_ranker.methods.base import fail_if_empty
from gene_ranker.dual_dataset import DualDataset

from pydeseq2.ds import DeseqStats, DeseqDataSet


@fail_if_empty
def deseq_shrinkage_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    case_cols = [x for x in dual_dataset.case.columns if x != dual_dataset.on]
    ctrl_cols = [x for x in dual_dataset.control.columns if x != dual_dataset.on]
    labels = ["case"] * len(case_cols) + ["control"] * len(ctrl_cols)
    metadata = pd.DataFrame({"sample": case_cols + ctrl_cols, "status": labels})
    metadata = metadata.set_index("sample")
    data = dual_dataset.merged.set_index(dual_dataset.on)

    data = round((2**data) - 1)
    data = data.applymap(int)
    data = data.transpose()

    data = DeseqDataSet(
        counts=data,
        metadata=metadata,
        design_factors="status",
        ref_level=["status", "control"],
        quiet=True,
    )

    data.deseq2()
    stats = DeseqStats(data)
    stats.summary()
    stats.lfc_shrink(coeff="status_case_vs_control")

    shrunk = stats.LFC

    result = pd.DataFrame(
        {dual_dataset.on: shrunk.index, "ranking": shrunk["status_case_vs_control"]}
    )

    return result
