from abc import ABC, abstractmethod
from Bio.SeqRecord import SeqRecord
from os import path
from glob import glob
from typing import Iterable, Optional, Sequence, Union

from .data import DnaSeq, ParseError

class SeqLoaderBase(ABC):

    @abstractmethod
    def load(self, file_path: str) -> Optional[Sequence[SeqRecord]]:
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

    def __load_file(self, file: str) -> Iterable[DnaSeq]:

        for loader in self.__loaders:
            results = loader.load(file)

            if results is not None:
                yield from (
                    DnaSeq(
                        seq=result,
                        seq_file=file
                    )
                    for result in results
                )

        # No loader exists for the given file
        return None

    def load_files(self) -> Iterable[SeqLoaderItem]:

        for filepath in self.__paths:
            pattern = path.join(filepath, "**", "*")

            for file in glob(pattern):
                try:
                    yield from self.__load_file(file)
                except ParseError as e:
                    yield e