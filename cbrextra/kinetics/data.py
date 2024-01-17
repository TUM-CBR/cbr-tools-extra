from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Tuple

class PartialInhibitionParameters(BaseModel):
    ksi: float
    km: float
    vmax: float
    beta: float

class PartialInhibitionRange(BaseModel):
    """Describes the range of values that can be assigned to
    the PartialInhibitionParameters. This gets interpreted as tensorflow
    regularizers that introduce a loss if the range is out of scope."""

    ksi: Tuple[Optional[float], Optional[float]] = (None, None)
    km: Tuple[Optional[float], Optional[float]] = (None, None)
    vmax: Tuple[Optional[float], Optional[float]] = (None, None)
    beta: Tuple[Optional[float], Optional[float]] = (None, None)

class ModelSpec(BaseModel):
    model_name : str
    model_parameters : PartialInhibitionParameters

class EvalArgs(BaseModel):
    model : ModelSpec
    data : List[float]

class Point2d(BaseModel):
    x: float
    y: float

class FitArgs(BaseModel):
    model : ModelSpec
    data : List[Point2d]
    iterations : int
    fit_range: Optional[PartialInhibitionRange] = None

class FitResult(BaseModel):
    model : ModelSpec
    original : ModelSpec

class EvalResult(BaseModel):
    results: List[Point2d]

class SimulationSpec(BaseModel):
    interval : int
    periods : int
    reps_per_unit : int

class SimulateArgs(BaseModel):
    model : ModelSpec
    simulation_spec: SimulationSpec
    initial_concentrations : List[float]

class FitSimulationArgs(BaseModel):
    model: ModelSpec
    simulation_spec: SimulationSpec
    data: Dict[float, List[float]]
    iterations : int
    fit_range: Optional[PartialInhibitionRange] = None

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