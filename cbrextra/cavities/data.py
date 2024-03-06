from pydantic import BaseModel, Field
from typing import List, Optional

class Points(BaseModel):
    id: str
    points: List[List[float]]

class FindCavitiesArgs(BaseModel):
    points_id: str
    min_volume: int
    max_volume: int

class InteractiveInput(BaseModel):
    find_cavities: Optional[FindCavitiesArgs] = Field(default=None)

class InteractiveOutput(BaseModel):
    pass