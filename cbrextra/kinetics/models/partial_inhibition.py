import keras

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