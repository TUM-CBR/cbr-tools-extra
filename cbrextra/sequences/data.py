from Bio.SeqRecord import SeqRecord
from typing import NamedTuple, Optional

class BlastEnv(NamedTuple):
    makeblastdb: str

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