from abc import abstractmethod, ABC
from enum import Enum
from typing import Dict, NamedTuple, Self

class Context(NamedTuple):
    cmdargs : dict

class ResultStatus(Enum):
    SUCCESS = 0
    NOT_REQUESTED = 1
    INVALID_ARGUMENTS = 2

class Result(NamedTuple):
    status :  ResultStatus

    @property
    def is_not_requested(self : Self) -> bool:
        return self.status == ResultStatus.NOT_REQUESTED

    @property
    def is_success(self : Self) -> bool:
        return self.status == ResultStatus.SUCCESS

    @property
    def is_requested(self : Self) -> bool:
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
    def success(self : Self) -> Result:
        return Result.success()

    @property
    def not_requested(self : Self) -> Result:
        return Result.not_requested()

    def invalid_arguments(self : Self, missing: Dict[str, str]) -> Result:
        return Result.invalid_arguments(missing)

    @abstractmethod
    def main(self, context : Context) -> Result:
        pass

class Operation(ABC):

    def __init__(self : Self, context : Context):
        self.__context = context

    @property
    def context(self : Self) -> Context:
        return self.__context