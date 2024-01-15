import tensorflow as tf
import keras

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

    def build(self, input_shape):

        self.ksi = self.add_weight(
            shape=(),
            initializer=constant_init(self.__ksi_0),
            trainable=True
        )

        self.beta = self.add_weight(
            shape=(),
            initializer=constant_init(self.__beta_0),
            trainable=True
        )

        self.km = self.add_weight(
            shape=(),
            initializer=constant_init(self.__km_0),
            trainable=True
        )

        self.vmax = self.add_weight(
            shape=(),
            initializer=constant_init(self.__vmax_0),
            trainable=True
        )

        self.one = self.add_weight(
            shape=(),
            initializer=constant_init(1),
            trainable=False
        )

    def call(self, inputs, *args, **kwargs):

        num = (self.one + self.beta * inputs / self.ksi) * inputs * self.vmax
        den = self.km + (1 + inputs / self.ksi) * inputs

        return num / den

    def to_model_spec(self, name: str) -> ModelSpec:
        return ModelSpec(
            model_name=name,
            model_parameters={
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
        self.__velocity_layer = PartialInhibitionLayer(ksi, km, vmax, beta)

    def call(self, inputs, training=None, mask=None):
        return self.__velocity_layer(inputs)

    def to_model_spec(self) -> ModelSpec:
        return self.__velocity_layer.to_model_spec(self.MODEL_NAME)

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
        self.__velocity_layer = PartialInhibitionLayer(ksi, km, vmax, beta)
        self.__simulation_spec = simulation_spec

    def to_model_spec(self) -> ModelSpec:
        return self.__velocity_layer.to_model_spec(self.MODEL_NAME)

    def call(self, inputs, training=None, mask=None):

        reps_per_unit = 10
        interval = self.__simulation_spec.interval
        periods = self.__simulation_spec.periods
        differential = tf.constant(1/(reps_per_unit * interval))
        result = []

        for i in range(0, periods):
            for j in range(0, reps_per_unit*interval):
                inputs -= differential * self.__velocity_layer(inputs)
            result.append(inputs)

        return tf.transpose(tf.convert_to_tensor(result))
