import tensorflow as tf
import keras
from keras import initializers

def constant_init(value: float):
    return {
        "class_name": "Constant",
        "config": {"value": value}
    }

class PartialInhibitionModel(keras.layers.Layer):

    MODEL_NAME = "partial_inhibition"

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

        dims = 1
        self.ksi = self.add_weight(
            shape=(dims,),
            initializer=constant_init(1/self.__ksi_0),
            trainable=True
        )

        self.beta = self.add_weight(
            shape=(dims,),
            initializer=constant_init(self.__beta_0),
            trainable=True
        )

        self.km = self.add_weight(
            shape=(dims,),
            initializer=constant_init(self.__km_0),
            trainable=True
        )

        self.vmax = self.add_weight(
            shape=(dims,),
            initializer=constant_init(self.__vmax_0),
            trainable=True
        )

        self.one = self.add_weight(
            shape=(dims,),
            initializer=constant_init(1),
            trainable=False
        )

    def call(self, inputs):

        num_1 = tf.matmul(
                tf.matmul(inputs, self.ksi),
                self.beta
            )

        num_2 = self.one = num_1

        num = self.matmul(
            num_2,
            input
        )

        num = (self.one + inputs * self.ksi * self.beta) * inputs

        den = self.km + \
            tf.matmul(
                1 + tf.matmul(inputs, self.ksi),
                input
            )

        return tf.matmul(num, tf.linalg.inv(den))