import pandas as pd
from typing import Optional


def two_way_in(x, y):
    """Check if all items in x are in y and vice-versa."""
    return all([k in y for k in x]) and all([k in x for k in y])


class DualDataset:
    """Represents two tangled datasets that can be dynamically updated

    The `case` and `control` DataFrames can always be found merged in the
    `merged` slot.
    If anything is updated (e.g. filtering), both case/control and merged
    slots are updated accordingly.

    The merge strategy is always `inner`.
    """

    def __init__(self, case: pd.DataFrame, control: pd.DataFrame, on="gene_id"):
        """Make a new DualDataset

        Args:
            case (pd.DataFrame): The case dataframe
            control (pd.DataFrame): The control dataframe
            on (str): The name of the column to merge on

        Raise:
            ValueError: If the 'on' column is not present in both case and control
                datasets.
            ValueError: If the two case/control datasets share more columns other
                than the 'on' column.
        """
        self.on: str = on

        if self.on not in case.columns or self.on not in control.columns:
            raise ValueError(f"Shared column '{self.on}' not in both datasets.")

        k = case.columns.tolist()
        k.remove(self.on)
        if any([x in control.columns for x in k]):
            raise ValueError("Case and control frames share columns other than `on`.")

        self._case: Optional[pd.DataFrame] = case.sort_values(by=self.on)
        self._case_cols = case.columns
        self._control: Optional[pd.DataFrame] = control.sort_values(by=self.on)
        self._control_cols = control.columns
        self._merged: Optional[pd.DataFrame] = None

    def sync(self):
        """Tangles the case/control datasets and then de-tangles them.

        This assures that the row order and quantity is the same among the two
        datasets as if they were just merged.
        """
        _ = self.merged
        self._case = None
        self._control = None

    def _detangle(self, columns) -> pd.DataFrame:
        cols = columns.tolist()
        return self.merged[cols]

    @property
    def merged(self):
        if self._merged is None:
            if self._case is None or self._control is None:
                raise ValueError("Cannot access merged without case and control")
            self._merged = self.case.merge(self.control, on=self.on)

        return self._merged

    @merged.setter
    def merged(self, value: pd.DataFrame):
        if not two_way_in(self._merged.columns, self.value.columns):
            raise ValueError("Cannot set new merged dataframe with different columns.")
        self._merged = value.sort_values(by=self.on)
        self._case = None
        self._control = None

    @property
    def case(self):
        if self._case is None:
            if self._merged is None:
                raise ValueError("Cannot access case without merged frame")
            self._case = self._detangle(self._case_cols)

        return self._case

    @property
    def control(self):
        if self._control is None:
            if self._merged is None:
                raise ValueError("Cannot access control without merged frame")
            self._control = self._detangle(self._control_cols)

        return self._control

    @case.setter
    def case(self, value):
        self._merged = None
        self._case_cols = value.columns
        self._case = value.sort_values(by=self.on)

    @control.setter
    def control(self, value):
        self._merged = None
        self._control_cols = value.columns
        self._control = value.sort_values(by=self.on)
