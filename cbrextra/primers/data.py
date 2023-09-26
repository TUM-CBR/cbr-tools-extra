from Bio.Seq import Seq
from pydantic import BaseModel
from typing import Dict, Iterable, List, Literal, NamedTuple

PrimerOrganism = Literal['E_COLI', 'P_PASTORIS']

class PrimerResult(NamedTuple):
    left_primer : str
    tm_left : float
    right_primer : str
    tm_right : float
    inner_seq : str
    tm_all : float

    def __get_penalty(self, target_tm : float, primer_tm : float):
        return abs(target_tm - primer_tm)

    def __get_error(
        self,
        target_tm : float,
        w_target_tm : float,
        primers_tm : float,
        w_primers_tm : float,
        w_tm_delta : float
    ) -> float:

        return sum([
            w_primers_tm*self.__get_penalty(primers_tm, self.tm_left),
            w_primers_tm*self.__get_penalty(primers_tm, self.tm_right),
            w_tm_delta*abs(self.tm_left - self.tm_right),
            w_target_tm*(self.tm_all and abs(self.tm_all - target_tm) or 0)
        ])

    @staticmethod
    def choose_best(
        opt1 : 'PrimerResult',
        opt2 : 'PrimerResult',
        target_tm : float,
        w_target_tm : float,
        primers_tm : float,
        w_primers_tm : float,
        w_tm_delta : float
    ) -> 'PrimerResult':
        score_opt1 = opt1.__get_error(
            target_tm,
            w_target_tm,
            primers_tm,
            w_primers_tm,
            w_tm_delta
        )
        score_opt2 = opt2.__get_error(
            target_tm,
            w_target_tm,
            primers_tm,
            w_primers_tm,
            w_tm_delta
        )
        if score_opt2 < score_opt1:
            return opt2
        else:
            return opt1

    @property
    def c_left_primer(self):
        return str(Seq(self.left_primer).complement())

    @property
    def c_right_primer(self):
        return str(Seq(self.right_primer).complement())

    @property
    def c_inner(self):
        return str(Seq(self.inner_seq).complement())

DesignPrimersResult = Dict[str, List[PrimerResult]]
DesignPrimersResults = Iterable[DesignPrimersResult]

class DesignPrimersArgs(BaseModel):
    sequence : str
    start : int
    codon_count : int
    min_length : int
    max_length : int
    organism : PrimerOrganism
