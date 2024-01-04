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

class FitArgs(BaseModel):
    model : ModelSpec
    data : List[Point2d]

class FitResult(BaseModel):
    model : ModelSpec
    original : ModelSpec

class EvalResult(BaseModel):
    results: List[Point2d]

class InteractiveInput(BaseModel):
    eval_model : Optional[EvalArgs] = Field(default=None)
    fit_model : Optional[FitArgs] = Field(default=None)

class InteractiveOutput(BaseModel):
    eval_result : Optional[EvalResult] = Field(default=None)
    fit_result : Optional[FitResult] = Field(default=None)