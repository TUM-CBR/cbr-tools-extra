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

    fit_range = args.fit_range

    if fit_range is not None:
        model.velocity_layer.with_range_loss(fit_range)

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
    model.fit(x_train, y_train, epochs=args.iterations)

    return FitResult(
        model = model.to_model_spec(),
        original = args.model
    )

def simulate_model(args: SimulateArgs) -> SimulateResult:
    model = select_simulation_model(args.model, args.simulation_spec)
    input = tf.convert_to_tensor(args.initial_concentrations)
    result = model(input)

    return SimulateResult(
        simulation_spec = args.simulation_spec,
        results = dict(
            (conc, list(map(float, values)))
            for conc, values in zip(args.initial_concentrations, result)
        )
    )

def fit_by_simulation(args: FitSimulationArgs):
    model = select_simulation_model(args.model, args.simulation_spec)

    fit_range = args.fit_range
    if fit_range is not None:
        model.velocity_layer.with_range_loss(fit_range)

    data = args.data
    x_train = tf.convert_to_tensor(list(data.keys()))
    y_train = tf.convert_to_tensor(list(data.values()))
    loss = MeanSquaredError()
    model.compile(
        optimizer='adam',
        loss=loss,
        metrics=['accuracy']
    )

    model.fit(x_train, y_train, epochs=args.iterations)

    return FitResult(
        model = model.to_model_spec(),
        original = args.model
    )