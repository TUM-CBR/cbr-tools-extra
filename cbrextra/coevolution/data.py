from Bio.Data.IUPACData import protein_letters
from pydantic import BaseModel
from typing import Dict, List, Optional, Sequence

K_OCCURRENCE_SCORE = 'occurrence'
K_EXLUSIVITY_SCORE = 'exclusivity'
K_CONFIDENCE_SCORE = 'confidence'
K_SYMMETRY_SCORE = 'symmetry'

TOTAL_RESIDUES = len(protein_letters)
RESIDUE_TO_INT = {
    res: i
    for i,res in enumerate(protein_letters)
}

IGAP = -1
RESIDUE_TO_INT['-'] = IGAP

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
    included_residues: Optional[Dict[int, List[str]]] = None

    def get_included(self, position: int) -> Optional[Sequence[int]]:

        included = self.included_residues
        if included is None or position not in included:
            return None
        
        residues = included[position]

        if len(residues) == 0:
            raise ValueError("At least one residue must be included per position.")

        def residue_to_number(residue: str) -> int:
            residue_number = RESIDUE_TO_INT.get(residue.upper())
            if residue_number is None:
                raise ValueError(f"The value '{residue}' is not a known residue.")
            
            return residue_number

        return [
            residue_to_number(residue)
            for residue in residues
        ]

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