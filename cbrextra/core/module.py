from abc import abstractmethod, ABC
from enum import Enum
from typing import Dict, NamedTuple

class Context(NamedTuple):
    cmdargs : dict

class ResultStatus(Enum):
    SUCCESS = 0
    NOT_REQUESTED = 1
    INVALID_ARGUMENTS = 2

class Result(NamedTuple):
    status :  ResultStatus

    @property
    def is_not_requested(self) -> bool:
        return self.status == ResultStatus.NOT_REQUESTED

    @property
    def is_success(self) -> bool:
        return self.status == ResultStatus.SUCCESS

    @property
    def is_requested(self) -> bool:
        return self.is_success

    @staticmethod
    def success():
        return Result(status=ResultStatus.SUCCESS)

    @staticmethod
    def not_requested():
        return Result(status=ResultStatus.NOT_REQUESTED)

    @staticmethod
    def invalid_arguments(missing : Dict[str, str]):
        return Result(status=ResultStatus.INVALID_ARGUMENTS)

class Module(ABC):

    @property
    def success(self) -> Result:
        return Result.success()

    @property
    def not_requested(self) -> Result:
        return Result.not_requested()

    def invalid_arguments(self, missing: Dict[str, str]) -> Result:
        return Result.invalid_arguments(missing)

    @abstractmethod
    def main(self, context : Context) -> Result:
        pass

class Operation(ABC):

    def __init__(self, context : Context):
        self.__context = context

    @property
    def context(self) -> Context:
        return self.__context