import shutil
import subprocess
import tempfile

import pandas as pd

from gene_ranker.dual_dataset import DualDataset
from gene_ranker.methods.base import (
    MissingExternalDependency,
    fail_if_empty,
    move_col_to_front,
)


@fail_if_empty
def cohen_d_ranking(dual_dataset: DualDataset) -> pd.DataFrame:
    if not shutil.which("fast-cohen"):
        raise MissingExternalDependency(
            (
                "Missing the fast-cohen executable."
                " See https://github.com/mrhedmad/fast-cohen/ to download."
            )
        )

    with (
        tempfile.NamedTemporaryFile() as case,
        tempfile.NamedTemporaryFile() as control,
        tempfile.NamedTemporaryFile() as result,
    ):
        move_col_to_front(dual_dataset.case, dual_dataset.on).to_csv(case, index=False)
        move_col_to_front(dual_dataset.control, dual_dataset.on).to_csv(
            control, index=False
        )

        subprocess.run(["fast-cohen", case.name, control.name, result.name])

        data = pd.read_csv(result)
        # Rename the 'row_names' col
        data = data.rename(columns={"row_names": dual_dataset.on, "cohen_d": "ranking"})

    return data
