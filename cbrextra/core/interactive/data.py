from enum import Enum
from pydantic import BaseModel, Field
from typing import List, Optional

class EntityType(Enum):
    STOP = 'stop'
    MESSAGE = 'message'

class InteractiveInput(BaseModel):
    uid : int
    entity_type: str
    payload: Optional[dict] = Field(default=None)

class ErrorCodes(Enum):
    UnknownError = 0
    JsonDecodeError = 1
    PydanticValidationError = 2
    UnknownInputType = 3

class InteractiveError(BaseModel):
    input_uids: List[int]
    error_code: ErrorCodes
    payload: Optional[str] = Field(default=None)
    message: str

class InteractiveValue(BaseModel):
    input_uids: List[int]
    payload: dict

class InteractiveOutput(BaseModel):
    uid : int
    value: Optional[InteractiveValue] = Field(default=None)
    error: Optional[InteractiveError] = Field(default=None)