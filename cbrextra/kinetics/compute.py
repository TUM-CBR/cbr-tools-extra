from keras.losses import MeanSquaredError
import tensorflow as tf
from typing import Any

from .data import *
from .models.factory import select_eval_model, select_simulation_model

def eval_model(args: EvalArgs) -> EvalResult:

    model = select_eval_model(args.model)
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
    model = select_eval_model(args.model)
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

def simulate_model(args: SimulateArgs) -> SimulateResult:
    model = select_simulation_model(args.model, args.simulation_spec)
    input = tf.convert_to_tensor(args.initial_concentrations)
    result = model(input)

    return SimulateResult(
        results = dict(
            (conc, list(result))
            for conc, values in zip(args.initial_concentrations, args.initial_concentrations)
        )
    )