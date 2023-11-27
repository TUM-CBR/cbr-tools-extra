from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from enum import Enum
from pydantic import BaseModel
from typing import cast, Dict, List, NamedTuple, Optional

class CascadeSequence(NamedTuple):
    seq_id : str
    seq : str

    @staticmethod
    def from_seq(seq: SeqRecord) -> 'CascadeSequence':

        if seq.id is None:
            raise ValueError("The sequence must have a valid id.")

        return CascadeSequence(
            seq_id = seq.id,
            seq = str(seq.seq)
        )


class CascadeStep(NamedTuple):
    step_id : int
    step_name : str
    sequences : List[CascadeSequence] = []

    @property
    def fasta(self) -> List[str]:
        return [
            f">{sequence.seq_id}\n{sequence.seq}"
            for sequence in self.sequences
        ]

    @staticmethod
    def from_fasta(file_location: str):
        return CascadeStep(
            step_id = 0,
            step_name=file_location,
            sequences= [
                CascadeSequence.from_seq(s)
                for s in SeqIO.parse(file_location, format='fasta')
            ]
        )


class Organism(NamedTuple):
    tax_id : int
    name: str

class CascadeStepOrganism(NamedTuple):
    organism : Organism
    identity : float
    sequence_match : str

default_include = [
    "Bacteria (taxid:2)",
    #"Archaea (taxid:2157)"
]
default_exclude = ["synthetic constructs (taxid:32630)"]

class FindOrganismsArgs(NamedTuple):
    step : CascadeStep
    excluded_organisms: Optional[List[str]] = None
    included_organisms: Optional[List[str]] = None
    num_results : int = 1000

def update_if_better(orgs: Dict[int, CascadeStepOrganism], org : CascadeStepOrganism):

    key = org.organism.tax_id
    prev = orgs.get(key)

    if prev is None \
        or org.identity > prev.identity:
        orgs[key] = org

class CascadeStepResult(NamedTuple):
    step : CascadeStep
    organisms : List[CascadeStepOrganism]

class CascadeReesult(NamedTuple):
    steps : List[CascadeStepResult]

class QueryStepPolicy(Enum):
    keep = "keep"
    replace = "replace"
    any = "any"

    @staticmethod
    def read(value : str) -> 'QueryStepPolicy':

        options = [option.value for option in QueryStepPolicy]

        if value in options:
            return cast(QueryStepPolicy, value)

        options_str = ", ".join(options)
        raise ValueError(f"Unexpected value '{value}', must be one of '{options_str}'")

class QueryCascadeStep(NamedTuple):
    step_id : int
    policy : QueryStepPolicy

class QueryCascadeArgs(NamedTuple):
    steps : List[QueryCascadeStep]
    max_identity_treshold : float

class OrganismResultEntry(BaseModel):
    tax_id : int
    name : str

class QueryCascadeResultStepEntry(BaseModel):
    step_id : int
    step_name : str
    identity : float

class QueryCascadeResultOrganismEntry(BaseModel):
    organism : OrganismResultEntry
    steps : List[QueryCascadeResultStepEntry]

class QueryCascadeResult(BaseModel):
    organisms : List[QueryCascadeResultOrganismEntry]