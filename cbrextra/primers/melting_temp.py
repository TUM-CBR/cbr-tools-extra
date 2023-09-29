import primer3
from primer3.bindings import DEFAULT_P3_ARGS

from .data import Primer3Args

class MeltingTemp:

    def __init__(
        self,
        args : Primer3Args
    ):
        self.__args = args

    def oligo_tm(self, seq : str) -> float:
        args = self.__args
        return primer3.calc_tm(
            seq,
            mv_conc=args.mv_conc or DEFAULT_P3_ARGS.mv_conc,
            dv_conc=args.dv_conc or DEFAULT_P3_ARGS.dv_conc,
            dntp_conc=args.dntp_conc or DEFAULT_P3_ARGS.dntp_conc,
            dna_conc=args.dna_conc or DEFAULT_P3_ARGS.dna_conc,
            annealing_temp_c=args.annealing_temp_c or DEFAULT_P3_ARGS.annealing_temp_c,
            max_nn_length=args.max_nn_length or DEFAULT_P3_ARGS.max_nn_length,
            dmso_conc=args.dmso_conc or DEFAULT_P3_ARGS.dmso_conc,
            dmso_fact=args.dmso_fact or DEFAULT_P3_ARGS.dmso_fact,
            formamide_conc=args.formamide_conc or DEFAULT_P3_ARGS.formamide_conc
        )
