from abc import ABC, abstractmethod
from os import path
from glob import glob
from typing import Iterable, Optional, Sequence, Union

from .data import DnaSeq, SequenceLoadException, SeqEntry

DnaSeqResult = Union[DnaSeq, SequenceLoadException]
SeqEntryResult = Union[SeqEntry, SequenceLoadException]

class SeqLoaderBase(ABC):

    @abstractmethod
    def load(self, file_path: str) -> Optional[Sequence[SeqEntryResult]]:
        raise NotImplementedError()
    
SeqLoaderItem = Union[DnaSeq, SequenceLoadException]

def to_result(item: SeqEntryResult, file: str) -> DnaSeqResult:

    if isinstance(item, SequenceLoadException):
        return item
    else:
        return DnaSeq(
            seq=item.seq,
            seq_file=file,
            tax_id=item.tax_id
        )

class SeqLoaderManager:

    def __init__(
        self,
        paths: Sequence[str],
        loaders: Sequence[SeqLoaderBase]
    ) -> None:
        
        self.__paths = paths
        self.__loaders = loaders

    def __load_file(self, file: str) -> Iterable[DnaSeqResult]:

        for loader in self.__loaders:
            results = loader.load(file)

            if results is not None:
                yield from (
                    to_result(result, file)
                    for result in results
                )

        # No loader exists for the given file
        return None

    def load_files(self) -> Iterable[SeqLoaderItem]:

        for filepath in self.__paths:
            pattern = path.join(filepath, "**", "*")

            for file in glob(pattern, recursive=True):
                try:
                    yield from self.__load_file(file)
                except SequenceLoadException as e:
                    yield e