from keras.layers import Layer
from keras.losses import MeanSquaredError
import tensorflow as tf
from typing import Any

from .data import *
from .models import PartialInhibitionModel

def select_model(model: ModelSpec) -> Layer:
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

def select_eval_model(args: EvalArgs) -> Layer:
    return select_model(args.model)

def eval_model(args: EvalArgs) -> EvalResult:

    model = select_eval_model(args)
    values = tf.convert_to_tensor(args.data)
    result : Any = model(values)
    
    return EvalResult(
        results = [
            Point2d(
                x = input,
                y = float(result)
            )
            for (input, result) in zip(values, result)
        ]
    )

def fit_model(args: FitArgs) -> FitResult:
    model = select_model(args.model)
    loss = MeanSquaredError()
    model.compile(
        optimizer='adam',
        loss=loss,
        metrics=['accuracy']
    )
    x_train = tf.convert_to_tensor([
        point.x
        for point in args.data
    ])
    y_train = tf.convert_to_tensor([
        point.y
        for point in args.data
    ])
    model.fit(x_train, y_train, epochs=5000)

    return FitResult(
        model = model.to_model_spec(),
        original = args.model
    )



