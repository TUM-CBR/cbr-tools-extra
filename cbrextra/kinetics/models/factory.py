from typing import Any

from ..data import ModelSpec, SimulationSpec
from .partial_inhibition import *

def select_eval_model(model: ModelSpec) -> Any:
    name = model.model_name
    parameters = model.model_parameters
    if name == PartialInhibitionModel.MODEL_NAME:
        return PartialInhibitionModel(
            ksi = parameters.ksi,
            km = parameters.beta,
            vmax = parameters.vmax,
            beta = parameters.beta
        )

    raise ValueError(f"Unknown model name '{name}'")

def select_simulation_model(model: ModelSpec, args: SimulationSpec) -> Any:
    name = model.model_name
    parameters = model.model_parameters
    if name == PartialInhibitionSimulationModel.MODEL_NAME:
        return PartialInhibitionSimulationModel(
            ksi = parameters.ksi,
            km = parameters.km,
            vmax = parameters.vmax,
            beta = parameters.beta,
            simulation_spec = args
        )

    raise ValueError(f"Unknown model name '{name}'")