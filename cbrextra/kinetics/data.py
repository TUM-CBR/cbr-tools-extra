from pydantic import BaseModel, Field
from typing import Dict, List, Optional

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

class InteractiveInput(BaseModel):
    eval_model : Optional[EvalArgs] = Field(default=None)

class InteractiveOutput(BaseModel):
    eval_result : Optional[EvalResult] = Field(default=None)