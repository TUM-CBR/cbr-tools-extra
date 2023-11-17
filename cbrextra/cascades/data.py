from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os import path
from typing import Dict, List, NamedTuple

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
    name : str
    sequences : List[CascadeSequence]

    @property
    def fasta(self) -> List[str]:
        return [
            f">{sequence.seq_id}\n{sequence.seq}"
            for sequence in self.sequences
        ]

    @staticmethod
    def from_fasta(file_location: str):
        return CascadeStep(
            name = path.basename(file_location),
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