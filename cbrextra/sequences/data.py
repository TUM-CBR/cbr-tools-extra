from Bio.SeqRecord import SeqRecord
from pydantic import BaseModel
from typing import List, NamedTuple, Optional, Sequence

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

class SearchArg(BaseModel):
    """
    A class that represents a search for the gnome of an organism.

    Attributes
    ----------
    search_id
        An identifier for this query. Search results will include this
        identyfier in order to match them with the query in the case
        of performing asynchronous searches.
    tax_ids
        A list of taxonomy ids for the organisms to consider.
    names
        Fuzzy names that will be used to identify the organisms by
        name.
    accession
        The id of the accession where the organism's genome was
        published.
    """
    search_id: str
    tax_ids: Optional[List[int]] = None
    names: Optional[List[str]] = None
    accession: Optional[str] = None

class SearchArgs(BaseModel):
    searches: List[SearchArg]

class SearchResultRecord(BaseModel):
    accession_id: str

class SearchResultError(BaseModel):
    message: str

class SearchResult(BaseModel):
    search_id: str
    records: List[SearchResultRecord]
    errors: List[SearchResultError]
