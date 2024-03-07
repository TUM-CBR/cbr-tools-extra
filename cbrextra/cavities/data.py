from igraph import Graph
from numpy import float64
from numpy.typing import NDArray
import pandas as pd
from pydantic import BaseModel, Field
from typing import Dict, List, NamedTuple, Optional

K_BOX_ID = "box_id"
K_BOX_SIZE = "size"
K_BOX_DEPTH = "depth"
K_BOX_X = "x"
K_BOX_Y = "y"
K_BOX_Z = "z"
K_BOX_CX = "cx"
K_BOX_CY = "cy"
K_BOX_CZ = "cz"

class Points(BaseModel):
    id: str
    points: List[List[float]]
    radii: List[float]

class FindCavitiesArgs(BaseModel):
    points_id: str
    min_volume: int
    max_volume: int

class InteractiveInput(BaseModel):
    find_cavities: Optional[List[FindCavitiesArgs]] = Field(default=None)

class CavityModel(BaseModel):
    points: List[List[float]]
    radii: List[float]

class CavitiesResult(BaseModel):
    cavities: Dict[str, List[CavityModel]]

class InteractiveOutput(BaseModel):
    cavities_result: Optional[CavitiesResult] = Field(default=None)

class ProteinCavity(NamedTuple):
    boxes: pd.DataFrame

    def to_cavity_model(self):
        points: NDArray[float64] = self.boxes[[K_BOX_X, K_BOX_Y, K_BOX_Z]].to_numpy()
        radii: NDArray[float64] = self.boxes[K_BOX_SIZE].to_numpy()
        return CavityModel(
            points = points.tolist(),
            radii = radii.tolist()
        )

class FindCavitiesGraph(NamedTuple):
    graph: Graph
    depth: int