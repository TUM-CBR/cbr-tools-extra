import primer3

class MeltingTemp:

    def __init__(
        self,
        dna_conc : float,
        dntp_conc : float
    ):

        self.__dna_conc = dna_conc
        self.__dntp_conc = dntp_conc

    def oligo_tm(self, seq : str) -> float:
        return primer3.calc_tm(
            seq,
            dna_conc=self.__dna_conc,
            dntp_conc=self.__dntp_conc
        )
