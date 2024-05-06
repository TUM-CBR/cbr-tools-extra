from abc import ABC, abstractmethod
from Bio.SeqRecord import SeqRecord
from os import path
from glob import glob
from typing import Iterable, Optional, Sequence, Union

from .data import DnaSeq, ParseError

class SeqLoaderBase(ABC):

    @abstractmethod
    def load(self, file_path: str) -> Optional[SeqRecord]:
        raise NotImplementedError()
    
SeqLoaderItem = Union[DnaSeq, ParseError]

class SeqLoaderManager:

    def __init__(
        self,
        paths: Sequence[str],
        loaders: Sequence[SeqLoaderBase]
    ) -> None:
        
        self.__paths = paths
        self.__loaders = loaders

    def __load_file(self, file: str) -> Optional[DnaSeq]:

        for loader in self.__loaders:
            result = loader.load(file)

            if result is not None:
                return DnaSeq(
                    seq=result,
                    seq_file=file
                )

        # No loader exists for the given file
        return None

    def load_files(self) -> Iterable[SeqLoaderItem]:

        for filepath in self.__paths:
            pattern = path.join(filepath, "**", "*")

            for file in glob(pattern):
                try:
                    result = self.__load_file(file)
                    if result is not None:
                        yield result
                except ParseError as e:
                    yield e