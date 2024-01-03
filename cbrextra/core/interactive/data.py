from enum import Enum
from pydantic import BaseModel
from typing import List, Optional

class InteractiveInput(BaseModel):
    uid : int
    entity_type: str
    payload: Optional[dict]

class ErrorCodes(Enum):
    UnknownError = 0
    JsonDecodeError = 1
    PydanticValidationError = 2

class InteractiveError(BaseModel):
    command_uids: List[int]
    error_code: ErrorCodes
    payload: Optional[str]
    message: str

class InteractiveValue(BaseModel):
    command_uids: List[int]
    payload: dict

class InteractiveOutput(BaseModel):
    uid : int
    value: Optional[InteractiveValue]
    error: Optional[InteractiveError]