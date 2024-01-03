from pydantic import BaseModel
from typing import Dict, List

class ModelSpec(BaseModel):
    model_name : str
    model_parameters : Dict[str, float]

class EvalArgs(BaseModel):
    model : ModelSpec
    data : List[float]

class Point2d(BaseModel):
    x: float
    y: float

class EvalResult(BaseModel):
    results: List[Point2d]
