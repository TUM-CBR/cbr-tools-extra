from Bio.SeqRecord import SeqRecord
from pydantic import BaseModel
from typing import NamedTuple, Optional, Sequence

class BlastEnv(NamedTuple):
    makeblastdb: str
    tblastn: str

class SequenceLoadException(Exception):

    def __init__(
        self,
        file: str,
        message: str
    ) -> None:
        super().__init__(f"Error reading '{file}': {message}")

        self.__file = file
        self.__message = message

    @property
    def file(self) -> str:
        return self.__file
    
    @property
    def message(self) -> str:
        return self.__message

class SeqEntry(NamedTuple):
    seq: SeqRecord
    tax_id: Optional[int] = None

class DnaSeq(NamedTuple):
    seq : SeqRecord
    seq_file : str
    tax_id: Optional[int] = None

class SequencesContext(NamedTuple):
    blast_context: BlastEnv
    db_file: str

    @property
    def db_blast_name(self):
        return f"{self.db_file}.blast"
    
class QueryOrganism(BaseModel):
    tax_id: int
    name: str

class QueryResult(BaseModel):
    id: str
    query_id: str
    organism: Optional[QueryOrganism]
    identity: float
    query_sequence: str
    hit_sequence: str
    mid_line: str
    file_location: Optional[str]

class QueryResults(BaseModel):
    results: Sequence[QueryResult]

class ErrorResult(BaseModel):
    message: str
    offending_file: str

class ErrorResults(BaseModel):
    results: Sequence[ErrorResult]