from pydantic import BaseModel
from typing import Any, Dict, List, Literal, Optional, Union

Scoring = Union[Literal['occurence'], Literal['symmetry']]

class Query(BaseModel):
    positions: List[int]
    max_results: int
    scoring: Scoring
    scoring_args: Dict[Any, Any]

class InteractiveInput(BaseModel):
    query: Optional[Query]

class CoevolutionEntry(BaseModel):
    residue: str
    score: float

class CoevolutionPosition(BaseModel):
    position: int
    by_position: Dict[int, CoevolutionEntry]

class CoevolutionResults(BaseModel):
    positions: Dict[int, CoevolutionPosition]

class InteractiveOutput(BaseModel):
    coevolution: Optional[CoevolutionResults]