from Bio.SeqRecord import SeqRecord
from typing import NamedTuple, Optional

class ParseError(Exception):
    pass

class DnaSeq(NamedTuple):
    seq : SeqRecord
    seq_file : str
    tax_id: Optional[int] = None