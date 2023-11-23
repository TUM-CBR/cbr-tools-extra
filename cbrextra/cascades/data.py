from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from enum import Enum
from typing import Dict, List, NamedTuple, Optional

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

class QueryCascadeStep(NamedTuple):
    step_id : int
    policy : QueryStepPolicy

class QueryCascadeArgs(NamedTuple):
    steps : List[QueryCascadeStep]
    max_identity_treshold : float

class QueryCascadeResultStepEntry(NamedTuple):
    step_id : int
    identity : float

class QueryCascadeResultOrganismEntry(NamedTuple):
    organism : Organism
    steps : List[QueryCascadeResultStepEntry]

class QueryCascadeResult(NamedTuple):
    organisms : List[QueryCascadeResultOrganismEntry]