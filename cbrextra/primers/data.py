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
    mv_conc: Optional[float] = None
    dv_conc: Optional[float] = None
    dntp_conc : Optional[float] = None
    dna_conc : Optional[float] = None
    annealing_temp_c : Optional[float] = None
    max_nn_length : Optional[int]  = None
    dmso_conc: Optional[float] = None
    dmso_fact: Optional[float] = None
    formamide_conc: Optional[float] = None

    @staticmethod
    def default() -> 'Primer3Args':
        return Primer3Args()

class DesignPrimersArgs(BaseModel):
    sequence : str
    start : int
    codon_count : int
    min_length : int
    max_length : int
    organism : PrimerOrganism
    primer3Args: Optional[Primer3Args]
