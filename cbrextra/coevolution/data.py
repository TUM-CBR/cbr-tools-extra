from pydantic import BaseModel
from typing import Dict, List, Optional

K_OCCURRENCE_SCORE = 'occurrence'
K_EXLUSIVITY_SCORE = 'exclusivity'
K_CONFIDENCE_SCORE = 'confidence'
K_SYMMETRY_SCORE = 'symmetry'

class Scoring(BaseModel):
    """
    Wether to residues are co-evolving or not is decided based
    on a score which constitutes many criteria. It is possible
    to weight each criteria differently by changing the paramters
    provided here.
    """

    occurrence_weight: float = 1
    exclusivity_weight: float = 1
    symmetry_weight: float = 1
    confidence_weight: float = 1
    confidence_treshold: float = 0.025
    confidence_center: float = 0.5
    confidence_concaveness: float = 8

    @property
    def __combined(self):
        return sum([
            self.exclusivity_weight,
            self.confidence_weight,
            self.symmetry_weight,
            self.occurrence_weight
        ])
    
    @property
    def occurrence_weight_scaled(self) -> float:
        return self.occurrence_weight / self.__combined

    @property
    def exclusivity_weight_scaled(self) -> float:
        return self.exclusivity_weight / self.__combined
    
    @property
    def confidence_weight_scaled(self) -> float:
        return self.confidence_weight / self.__combined
    
    @property
    def symmetry_weight_scaled(self) -> float:
        return self.symmetry_weight / self.__combined

class Query(BaseModel):
    positions: List[int]
    max_results: int
    scoring: Scoring = Scoring()

class InteractiveRequest(BaseModel):
    query: Optional[Query]

class CoevolutionEntry(BaseModel):
    residue_1: str
    residue_2: str
    score: float
    score_occurence: float
    score_exclusivity: float
    score_confidence: float
    score_symmetry: float

class CoevolutionPosition(BaseModel):
    position: int
    by_position: Dict[int, CoevolutionEntry]

class CoevolutionResults(BaseModel):
    positions: Dict[int, CoevolutionPosition]

class InteractiveResponse(BaseModel):
    coevolution: Optional[CoevolutionResults]