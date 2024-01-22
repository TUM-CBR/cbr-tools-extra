import tensorflow as tf
import keras
from keras.constraints import NonNeg
from keras.regularizers import Regularizer
from typing import Any, Optional

from ..data import ModelSpec, PartialInhibitionParameters, PartialInhibitionRange, SimulationSpec

def constant_init(value: float):
    return {
        "class_name": "Constant",
        "config": {"value": value}
    }

class RangeRegularizer(Regularizer):

    def __init__(self, low: Optional[float], high: Optional[float]):
        self.low_in = low
        self.low = tf.constant(
            0 if low is None else low,
            dtype=tf.float32,
            shape=()
        )
        self.low_w = tf.constant(
            0 if low is None else 1,
            dtype=tf.float32,
            shape=()
        )

        self.high_in = high
        self.high = tf.constant(
            0 if high is None else high,
            dtype=tf.float32,
            shape=()
        )
        self.high_w = tf.constant(
            0 if high is None else 1,
            dtype=tf.float32,
            shape=()
        )

        self.zero = tf.constant(0, dtype=tf.float32, shape=())

    @tf.function
    def __call__(self, value: Any) -> Any:

        low_loss = tf.cond(
            pred = tf.less(value, self.low),
            true_fn = lambda: tf.exp(tf.abs(self.low - value)) * self.low_w,
            false_fn = lambda: self.zero
        )

        high_loss = tf.cond(
            pred = tf.greater(value, self.high),
            true_fn = lambda: tf.exp(tf.abs(value - self.high)) * self.high_w,
            false_fn = lambda: self.zero
        )

        return low_loss + high_loss
    
    def get_config(self):
        return {
            'low': self.low_in,
            'high': self.high_in
        }

class PartialInhibitionLayer(keras.layers.Layer):

    def __init__(
        self,
        ksi : float,
        km : float,
        vmax : float,
        beta : float,
        ranges: Optional[PartialInhibitionRange] = None
    ):
        super().__init__()

        self.__ksi_0 = ksi
        self.__km_0 = km
        self.__vmax_0 = vmax
        self.__beta_0 = beta
        self.ranges = ranges

    def build(self, input_shape : Any):

        ranges = self.ranges
        if ranges is not None:
            ksi_reg = RangeRegularizer(*ranges.ksi)
            beta_reg = RangeRegularizer(*ranges.beta)
            km_reg = RangeRegularizer(*ranges.km)
            vmax_reg = RangeRegularizer(*ranges.vmax)
        else:
            ksi_reg = None
            beta_reg = None
            km_reg = None
            vmax_reg = None

        self.ksi : Any = self.add_weight(
            shape=(),
            initializer=constant_init(self.__ksi_0),
            trainable=True,
            constraint=NonNeg(),
            regularizer=ksi_reg
        )
        self.beta : Any = self.add_weight(
            shape=(),
            initializer=constant_init(self.__beta_0),
            trainable=True,
            constraint=NonNeg(),
            regularizer=beta_reg
        )
#
        self.km : Any = self.add_weight(
            shape=(),
            initializer=constant_init(self.__km_0),
            trainable=True,
            constraint=NonNeg(),
            regularizer=km_reg
        )
#
        self.vmax : Any = self.add_weight(
            shape=(),
            initializer=constant_init(self.__vmax_0),
            trainable=True,
            constraint=NonNeg(),
            regularizer=vmax_reg
        )


    def call(self, inputs : Any, training = False, *args : Any, **kwargs : Any) -> Any:

        epsillon = tf.constant(0.00000000000019236959237, dtype=tf.float32, shape=())

        # adding epsillon to avoid a division by 0
        ksi = self.ksi + epsillon
        num = 1 + ((self.beta * inputs) / ksi)
        den = self.km + (inputs * ((1 + inputs) / ksi))

        return inputs * self.vmax * num / (den + epsillon)

    def to_model_spec(self, name: str) -> ModelSpec:
        return ModelSpec(
            model_name = name,
            model_parameters = PartialInhibitionParameters(
                ksi = float(self.ksi),
                km = float(self.km),
                vmax = float(self.vmax),
                beta = float(self.beta)
            )
        )

class PartialInhibitionModel(keras.Model):

    MODEL_NAME = "partial_inhibition"

    def __init__(
        self,
        ksi : float,
        km : float,
        vmax : float,
        beta : float,
        ranges: Optional[PartialInhibitionRange] = None
    ):
        super().__init__()
        self.velocity_layer = PartialInhibitionLayer(ksi, km, vmax, beta, ranges)

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
        simulation_spec: SimulationSpec,
        ranges: Optional[PartialInhibitionRange] = None
    ):
        super().__init__()
        self.velocity_layer = PartialInhibitionLayer(ksi, km, vmax, beta, ranges)
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
        self.reps_per_period = tf.constant(reps_per_unit * interval, dtype=tf.int32, shape=())
        
        self.differential = tf.constant(1/(reps_per_unit * interval), dtype=tf.float32)
        #self.initial_conc = tf.constant(0, dtype=tf.float32, shape=input_shape)

    @tf.function
    def simulation_step(
            self,
            p: Any,
            product_concentrations: Any,
            substrate_concentrations: Any,
            result: Any
        ) -> Any:

        (p, _1, p_conc, s_conc, _2, result) = tf.while_loop(
            lambda p, step, p_conc, s_conc, rate, result: tf.less(step, self.reps_per_period),
            lambda p, step, p_conc, s_conc, rate, result: (
                p,
                step +1,
                p_conc + rate,
                s_conc - rate,
                self.differential * self.velocity_layer(s_conc - rate), result
            ),
            (
                p,
                tf.constant(0, dtype=tf.int32, shape=()),
                product_concentrations,
                substrate_concentrations,
                self.differential * self.velocity_layer(substrate_concentrations),
                result
            ),
            parallel_iterations=1
        )

        return (p + 1, p_conc, s_conc, result.write(p, p_conc))

    def call(self, inputs : Any, training : Any = None, mask : Any = None) -> Any:
        return self.simulate(inputs)

    @tf.function
    def simulate(self, inputs: Any) -> Any:
        (_1,_2,_3,result) = tf.while_loop(
            lambda p,p_conc,s_conc,result: tf.less(p, self.periods),
            self.simulation_step,
            (
                tf.constant(0, dtype=tf.int32, shape=()),
                inputs*0,
                inputs,
                tf.TensorArray(tf.float32, size=self.periods, clear_after_read=False, dynamic_size=False)
            ),
            parallel_iterations=1
        )
        return tf.transpose(result.stack())