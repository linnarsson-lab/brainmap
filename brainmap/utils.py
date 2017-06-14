from typing import *
import numpy as np
from collections import OrderedDict


class LimitedSizeDict(OrderedDict):
    def __init__(self, *args: Any, **kwds: Any) -> None:
        self.size_limit = kwds.pop("size_limit", None)  # type: int
        OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key: Any, value: Any) -> None:
        OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self) -> None:
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)


def one_hot_encoding(array2d: np.ndarray) -> List[np.ndarray]:
    temp = np.array(array2d, dtype=int)
    labels, temp.flat[:] = np.unique(temp, return_inverse=True)
    arrays = [np.zeros(temp.shape, dtype=bool) for _ in range(len(labels))]
    for i in range(temp.shape[0]):
        for j in range(temp.shape[1]):
            arrays[temp[i, j]][i, j] = True
    return arrays
