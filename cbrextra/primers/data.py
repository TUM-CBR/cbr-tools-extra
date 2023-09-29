from Bio.Seq import Seq
from pydantic import BaseModel
from typing import Dict, Iterable, Literal, Optional

PrimerOrganism = Literal['E_COLI', 'P_PASTORIS']

class PrimerResult(BaseModel):
    left_primer : str
    tm_left : float
    right_primer : str
    tm_right : float
    inner_seq : str
    tm_all : float
    position : int
    amino_acid : str

    @property
    def c_left_primer(self):
        return str(Seq(self.left_primer).complement())

    @property
    def c_right_primer(self):
        return str(Seq(self.right_primer).complement())

    @property
    def c_inner(self):
        return str(Seq(self.inner_seq).complement())

class DesignPrimersResults(BaseModel):
    primers : Iterable[PrimerResult]
    plasmid : str
    codon_mappings : Dict[str, str]

class Primer3Args(BaseModel):
    mv_conc: Optional[float]
    dv_conc: Optional[float]
    dntp_conc : Optional[float]
    dna_conc : Optional[float]
    annealing_temp_c : Optional[float]
    max_nn_length : Optional[int]
    dmso_conc: Optional[float]
    dmso_fact: Optional[float]
    formamide_conc: Optional[float]

    @staticmethod
    def default() -> 'Primer3Args':
        return Primer3Args(
            mv_conc=None,
            dv_conc=None,
            dna_conc=None,
            dntp_conc=None,
            annealing_temp_c=None,
            max_nn_length=None,
            dmso_conc=None,
            dmso_fact=None,
            formamide_conc=None
        )

class DesignPrimersArgs(BaseModel):
    sequence : str
    start : int
    codon_count : int
    min_length : int
    max_length : int
    organism : PrimerOrganism
    primer3Args: Optional[Primer3Args]
