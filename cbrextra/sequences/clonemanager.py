from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from enum import Enum
from os import path
from typing import Optional

from .data import ParseError
from .sequence import SeqLoaderBase

DNA_BASES = [b'C', b'T', b'A', b'G']
MIN_LEGIT_LENGTH = 100
CM5_EXT = '.cm5'

class State(Enum):
    Skipping = 0
    Reading = 1
    Done = 3

class CMLoader(SeqLoaderBase):

    def load(self, file_path: str) -> Optional[SeqRecord]:

        if file_path[-4:].lower() != CM5_EXT:
            return None
        
        seq = self.read_seq(file_path)

        if seq is None:
            raise ParseError("No sequence found in file!")
        

        return SeqRecord(
            seq=Seq(seq),
            id=path.basename(file_path)[0:-4]
        )

    def read_seq(self, file_path: str) -> Optional[bytes]:
        
        with open(file_path, 'rb') as stream:

            seq = b''
            state = State.Skipping

            while state != State.Done:
                
                value = stream.read(1)

                if value in DNA_BASES and state == State.Reading:
                    seq += value
                elif value in DNA_BASES and state == State.Skipping:
                    seq += value
                    state = State.Reading
                elif value not in DNA_BASES \
                    and state == State.Reading \
                    and len(seq) >= MIN_LEGIT_LENGTH:
                    return seq
                else:
                    state = State.Skipping
                    seq = b''

                if len(value) == 0:
                    state = State.Done

        return None