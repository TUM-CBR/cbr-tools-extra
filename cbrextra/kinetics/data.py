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

class SimulationSpec(BaseModel):
    interval : int
    periods : int

class SimulateArgs(BaseModel):
    model : ModelSpec
    simulation_spec: SimulationSpec
    initial_concentrations : List[float]

class FitSimulationArgs(BaseModel):
    model: ModelSpec
    simulation_spec: SimulationSpec
    data: Dict[float, List[float]]

class SimulateResult(BaseModel):
    simulation_spec: SimulationSpec
    results: Dict[float, List[float]]

class InteractiveInput(BaseModel):
    eval_model : Optional[EvalArgs] = Field(default=None)
    fit_model : Optional[FitArgs] = Field(default=None)
    simulate_model : Optional[SimulateArgs] = Field(default=None)
    fit_simulation : Optional[FitSimulationArgs] = Field(default=None)

class InteractiveOutput(BaseModel):
    eval_result : Optional[EvalResult] = Field(default=None)
    fit_result : Optional[FitResult] = Field(default=None)
    simulate_result : Optional[SimulateResult] = Field(default = None)