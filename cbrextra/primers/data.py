from Bio.Seq import Seq
from pydantic import BaseModel
from typing import Dict, Iterable, Literal

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

class DesignPrimersArgs(BaseModel):
    sequence : str
    start : int
    codon_count : int
    min_length : int
    max_length : int
    organism : PrimerOrganism
