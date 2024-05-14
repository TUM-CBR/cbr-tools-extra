from abc import ABC, abstractmethod
from glob import glob
import json
from json import JSONDecodeError
from os import path
from typing import Any, Dict, Iterable, NamedTuple, Optional, Sequence, Tuple, Union

from .data import DnaSeq, SequenceLoadException, SeqEntry

DnaSeqResult = Union[DnaSeq, SequenceLoadException]
SeqEntryResult = Union[SeqEntry, SequenceLoadException]

CONFIG_EXT = "sco.json"

class ConfigEntry(NamedTuple):
    configuration: Dict[Any, Any]
    path: str

class ConfigContext(NamedTuple):
    configurations: Sequence[ConfigEntry]

class SeqLoaderBase(ABC):

    @abstractmethod
    def load(self, file_path: str) -> Optional[Sequence[SeqEntryResult]]:
        raise NotImplementedError()
    
class SeqLoaderFactoryBase(ABC):

    @abstractmethod
    def configure(self, config: ConfigContext) -> SeqLoaderBase:
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
        loaders: Sequence[SeqLoaderFactoryBase]
    ) -> None:
        
        self.__paths = paths
        self.__loaders = loaders

    def __load_file(self, loaders: Sequence[SeqLoaderBase], file: str) -> Iterable[DnaSeqResult]:

        for loader in loaders:
            results = loader.load(file)

            if results is not None:
                yield from (
                    to_result(result, file)
                    for result in results
                )

        # No loader exists for the given file
        return None
    
    def __load_configurations(self) -> Iterable[Union[SequenceLoadException, ConfigEntry]]:

        for filepath in self.__paths:
            pattern = path.join(filepath, "**", f"*.{CONFIG_EXT}")

            for file in glob(pattern,  recursive=True):

                try:
                    with open(file, 'r') as config_stream:
                        json_config = json.load(config_stream)
                        yield ConfigEntry(
                            configuration=json_config,
                            path=file
                        )
                except JSONDecodeError as e:
                    yield SequenceLoadException(
                        file,
                        str(e)
                    )

    def __configure_loaders(self) -> Tuple[Sequence[SeqLoaderBase], Sequence[SequenceLoadException]]:

        configurations = self.__load_configurations()
        configuration_context = ConfigContext(
            configurations = [
                config
                for config in configurations
                    if isinstance(config, ConfigEntry)
            ]
        )

        return (
            [loader.configure(configuration_context) for loader in self.__loaders],
            [exn for exn in configurations if isinstance(exn, SequenceLoadException)]
        )

    def load_files(self) -> Iterable[SeqLoaderItem]:

        (loaders, config_errors) = self.__configure_loaders()

        for error in config_errors:
            yield error

        for filepath in self.__paths:
            pattern = path.join(filepath, "**", "*")

            for file in glob(pattern, recursive=True):
                try:
                    yield from self.__load_file(loaders, file)
                except SequenceLoadException as e:
                    yield e