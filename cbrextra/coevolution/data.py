from pydantic import BaseModel
from typing import Any, Dict, List, Optional

K_OCCURRENCE_SCORE = 'occurrence'
K_EXLUSIVITY_SCORE = 'exclusivity'
K_CONSERVED_SCORE = 'conserved'

class Scoring(BaseModel):
    """
    Wether to residues are co-evolving or not is decided based
    on a score which constitutes many criteria. It is possible
    to weight each criteria differently by changing the paramters
    provided here.
    """

    occurence_weight: float = 1
    exclusivity_weight: float = 1
    conserved_weight: float = 1

class Query(BaseModel):
    positions: List[int]
    max_results: int
    scoring: Scoring = Scoring()

class InteractiveInput(BaseModel):
    query: Optional[Query]

class CoevolutionEntry(BaseModel):
    residue_1: str
    residue_2: str
    score: float
    score_occurence: float
    score_exclusivity: float
    score_conserved: float

class CoevolutionPosition(BaseModel):
    position: int
    by_position: Dict[int, CoevolutionEntry]

class CoevolutionResults(BaseModel):
    positions: Dict[int, CoevolutionPosition]

class InteractiveOutput(BaseModel):
    coevolution: Optional[CoevolutionResults]