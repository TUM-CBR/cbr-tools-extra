from numpy import float32
import tensorflow as tf
import keras
from typing import Any

from ..data import ModelSpec, SimulationSpec

def constant_init(value: float):
    return {
        "class_name": "Constant",
        "config": {"value": value}
    }

class PartialInhibitionLayer(keras.layers.Layer):

    def __init__(
        self,
        ksi : float,
        km : float,
        vmax : float,
        beta : float
    ):
        super().__init__()

        self.__ksi_0 = ksi
        self.__km_0 = km
        self.__vmax_0 = vmax
        self.__beta_0 = beta

    def build(self, input_shape : Any):

        self.ksi : Any = self.add_weight(
            shape=(),
            initializer=constant_init(self.__ksi_0),
            trainable=True
        )

        self.beta : Any = self.add_weight(
            shape=(),
            initializer=constant_init(self.__beta_0),
            trainable=True
        )

        self.km : Any = self.add_weight(
            shape=(),
            initializer=constant_init(self.__km_0),
            trainable=True
        )

        self.vmax : Any = self.add_weight(
            shape=(),
            initializer=constant_init(self.__vmax_0),
            trainable=True
        )

        self.one : Any = self.add_weight(
            shape=(),
            initializer=constant_init(1),
            trainable=False
        )

    def call(self, inputs : Any, *args : Any, **kwargs : Any) -> Any:

        num = (self.one + self.beta * inputs / self.ksi) * inputs * self.vmax
        den = self.km + (1 + inputs / self.ksi) * inputs

        return num / den

    def to_model_spec(self, name: str) -> ModelSpec:
        return ModelSpec(
            model_name = name,
            model_parameters = {
                "ksi": float(self.ksi),
                "km": float(self.km),
                "vmax": float(self.vmax),
                "beta": float(self.beta)
            }
        )

class PartialInhibitionModel(keras.Model):

    MODEL_NAME = "partial_inhibition"

    def __init__(
        self,
        ksi : float,
        km : float,
        vmax : float,
        beta : float
    ):
        super().__init__()
        self.velocity_layer = PartialInhibitionLayer(ksi, km, vmax, beta)

    def call(self, inputs : Any, training : Any = None, mask : Any = None) -> Any:
        return self.velocity_layer(inputs)

    def to_model_spec(self) -> ModelSpec:
        return self.velocity_layer.to_model_spec(self.MODEL_NAME)

class PartialInhibitionSimulationModel(keras.Model):

    MODEL_NAME = "partial_inhibition"

    def __init__(
        self,
        ksi : float,
        km : float,
        vmax : float,
        beta : float,
        simulation_spec: SimulationSpec
    ):
        super().__init__()
        self.velocity_layer = PartialInhibitionLayer(ksi, km, vmax, beta)
        self.__simulation_spec = simulation_spec

    def to_model_spec(self) -> ModelSpec:
        return self.velocity_layer.to_model_spec(self.MODEL_NAME)
    
    def build(self, input_shape : Any):
        self.periods = tf.constant(
            self.__simulation_spec.periods,
            dtype=tf.int32,
            shape=()
        )

        reps_per_unit = self.__simulation_spec.reps_per_unit
        interval = self.__simulation_spec.interval
        self.differential = tf.constant(1/(reps_per_unit * interval), dtype=tf.float32)
        self.reps_per_period = tf.constant(reps_per_unit * interval, dtype=tf.int32)
        self.initial_conc = tf.constant(0, dtype=tf.float32, shape=input_shape)

    def call(self, inputs : Any, training : Any = None, mask : Any = None):
        result = tf.TensorArray(tf.float32, size=self.periods, clear_after_read=False, dynamic_size=False)
        product_concentrations = self.initial_conc
        substrate_concentrations = inputs

        for p in tf.range(0, self.periods):

            for _ in tf.range(0, self.reps_per_period):
                rate = self.differential * self.velocity_layer(substrate_concentrations)
                substrate_concentrations -= rate
                product_concentrations += rate

            result = result.write(p, product_concentrations)

        return tf.transpose(result.stack())
