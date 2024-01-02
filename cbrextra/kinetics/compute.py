from keras.layers import Layer

from .data import *
from .models import PartialInhibitionModel

def select_eval_model(args: EvalArgs) -> Layer:

    model = args.model
    name = model.model_name
    parameters = model.model_parameters
    if name == PartialInhibitionModel.MODEL_NAME:
        return PartialInhibitionModel(
            ksi = parameters["ksi"],
            km = parameters["km"],
            vmax = parameters["vmax"],
            beta = parameters["beta"]
        )

    raise ValueError(f"Unknown model name '{name}'")

def eval_model(args: EvalArgs) -> EvalResult:

    model = select_eval_model(args)
    return model(args.data)
    
    return EvalResult(
        results = [
            Point2d(
                x = value,
                y = model([value])
            )
            for value in args.data
        ]
    )