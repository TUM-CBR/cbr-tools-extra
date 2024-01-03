from abc import ABC, abstractmethod
from typing import Generic, NamedTuple, TypeVar

TMessageIn = TypeVar('TMessageIn')
TMessageOut = TypeVar('TMessageOut')

class MessageContetBase(ABC, Generic[TMessageIn, TMessageOut]):
    pass

class ParseMessageArgs(NamedTuple, Generic[TMessageIn, TMessageOut]):
    pass

class SerializeMessageArgs(NamedTuple, Generic[TMessageIn, TMessageOut]):
    pass

class OnMessageContextBase(MessageContetBase[TMessageIn, TMessageOut]):

    @abstractmethod
    def send(
        self,
        message: TMessageOut
    ) -> None:
        raise NotImplementedError()

class OnMessageArgs(NamedTuple, Generic[TMessageIn, TMessageOut]):
    context: OnMessageContextBase[TMessageIn, TMessageOut]
    message: TMessageIn

class InteractiveSpec(ABC, Generic[TMessageIn, TMessageOut]):

    @abstractmethod
    def parse_message(
        self,
        args: 'ParseMessageArgs[TMessageIn, TMessageOut]'
    ) -> TMessageIn:
        raise NotImplementedError()

    @abstractmethod
    def serialize_message(
        self,
        args: 'SerializeMessageArgs[TMessageIn, TMessageOut]'
    ) -> dict:
        raise NotImplementedError()

    def on_message(
        self,
        args: 'OnMessageArgs[TMessageIn, TMessageOut]'
    ):
        raise NotImplementedError()