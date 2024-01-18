import tensorflow as tf
import keras
from typing import Any, Union

from ..data import ModelSpec, PartialInhibitionParameters, PartialInhibitionRange, SimulationSpec

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
            model_parameters = PartialInhibitionParameters(
                ksi = float(self.ksi),
                km = float(self.km),
                vmax = float(self.vmax),
                beta = float(self.beta)
            )
        )

    @tf.function
    def divergence_penalty(self, upper: Any, lower: Any) -> Any:
        """Penalty applied when a variable goes outside a desried range.
        The difference between the current value and the boundry is multiplied
        by a distortion factor to amplify the effect. The resulting value is
        raised to a power of 4. This has been chosen as the SDG algorithm will
        look at the first two derivatives to guide its path. A power of 4 is
        the smallest power where the 1st and 2nd derivative will not be constant."""

        return tf.pow(
            tf.multiply(
                self.loss_distortion,
                tf.abs(upper - lower)
            ),
            self.divergence_penalty_power
        )

    def min_value_loss(
        self,
        variable: Any,
        min_value: Any
    ) -> Any:
        return tf.cond(
            pred = tf.greater_equal(variable, min_value),
            true_fn = tf.function(lambda: self.no_penalty),
            false_fn = tf.function(lambda: self.divergence_penalty(min_value, variable))
        )

    @tf.function
    def max_value_loss(
        self,
        variable: Any,
        max_value: Any
    ) -> Any:

        return tf.cond(
            pred = tf.less_equal(variable, max_value),
            true_fn = tf.function(lambda: self.no_penalty),
            false_fn = tf.function(lambda: self.divergence_penalty(variable, max_value))
        )
    
    def with_range_loss(
        self,
        ranges: PartialInhibitionRange
    ):
        self.loss_distortion = tf.constant(100000, dtype=tf.float32, shape=())
        self.divergence_penalty_power = tf.constant(4, dtype=tf.float32, shape=())
        self.no_penalty = tf.constant(0, dtype=tf.float32, shape=())

        min_ksi, max_ksi = ranges.ksi
        if min_ksi is not None:
            self.min_ksi = tf.constant(min_ksi, dtype=tf.float32, shape=())
            self.add_loss(lambda: self.min_value_loss(self.ksi, self.min_ksi))

        if max_ksi is not None:
            self.max_ksi = tf.constant(max_ksi, dtype=tf.float32, shape=())
            self.add_loss(lambda: self.max_value_loss(self.ksi, self.max_ksi))

        min_km, max_km = ranges.km
        if min_km is not None:
            self.min_km = tf.constant(min_km, dtype=tf.float32, shape=())
            self.add_loss(lambda: self.min_value_loss(self.km, self.min_km))

        if max_km is not None:
            self.max_km = tf.constant(max_km, dtype=tf.float32, shape=())
            self.add_loss(lambda: self.max_value_loss(self.km, self.max_km))

        min_vmax, max_vmax = ranges.vmax
        if min_vmax is not None:
            self.min_vmax = tf.constant(min_vmax, dtype=tf.float32, shape=())
            self.add_loss(lambda: self.min_value_loss(self.vmax, self.min_vmax))

        if max_vmax is not None:
            self.max_vmax = tf.constant(max_vmax, dtype=tf.float32, shape=())
            self.add_loss(lambda: self.max_value_loss(self.vmax, self.max_vmax))

        min_beta, max_beta = ranges.beta
        if min_beta is not None:
            self.min_beta = tf.constant(min_beta, dtype=tf.float32, shape=())
            self.add_loss(lambda: self.min_value_loss(self.beta, self.min_beta))

        if max_beta is not None:
            self.max_beta = tf.constant(max_beta, dtype=tf.float32, shape=())
            self.add_loss(lambda: self.max_value_loss(self.beta, self.max_beta))

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

    @tf.function
    def call(self, inputs : Any, training : Any = None, mask : Any = None) -> Any:
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