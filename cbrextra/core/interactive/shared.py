from abc import ABC, abstractmethod
from typing import Any, Dict, Generic, List, NamedTuple, Optional, TypeVar

TMessageIn = TypeVar('TMessageIn')
TMessageOut = TypeVar('TMessageOut')

class InputHeader(NamedTuple):
    uid : int

class MessageContextBase(ABC, Generic[TMessageIn, TMessageOut]):
    pass

class ParseMessageArgs(NamedTuple, Generic[TMessageIn, TMessageOut]):
    payload : Dict[Any, Any]

class SerializeMessageArgs(NamedTuple, Generic[TMessageIn, TMessageOut]):
    value : TMessageOut

class OnMessageContextBase(MessageContextBase[TMessageIn, TMessageOut]):

    @abstractmethod
    def send(
        self,
        message: TMessageOut,
        message_ids : List[int]
    ) -> None:
        raise NotImplementedError()

class OnMessageArgs(NamedTuple, Generic[TMessageIn, TMessageOut]):
    context: OnMessageContextBase[TMessageIn, TMessageOut]
    header: InputHeader 
    message: TMessageIn

    def send(
        self,
        message : TMessageOut,
        message_ids : Optional[List[int]] = None
    ):
        message_ids = \
            [self.header.uid] if message_ids is None \
            else message_ids
        self.context.send(
            message,
            message_ids
        )

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
    ) -> Dict[Any, Any]:
        raise NotImplementedError()

    @abstractmethod
    def on_message(
        self,
        args: 'OnMessageArgs[TMessageIn, TMessageOut]'
    ) -> None:
        raise NotImplementedError()