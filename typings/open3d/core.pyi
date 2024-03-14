from enum import Enum
from numpy import float64
from numpy.typing import ArrayLike, NDArray

class Dtype(Enum):
    Float32 = ...

class Tensor:

    def __init__(self, values: ArrayLike, dtype: Dtype = ...) -> None: ...

    def numpy(self) -> NDArray[float64]: ...