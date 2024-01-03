import json
from pydantic import ValidationError
from typing import TextIO

from .data import *
from .shared import *

class InteractiveHandler(Generic[TMessageIn, TMessageOut]):

    class OnMessageContext(OnMessageContextBase):

        def __init__(
            self,
            handler: 'InteractiveHandler'
        ) -> None:
            super().__init__()
            self.__handler = handler

        def send(
            self,
            message: TMessageOut
        ) -> None:
            pass

    def __init__(
        self,
        spec: InteractiveSpec[TMessageIn, TMessageOut],
        input_stream: TextIO,
        output_stream: TextIO
    ):
        self.__spec = spec
        self.__input_stream = input_stream
        self.__output_stream = output_stream
        self.__on_message_context = self.OnMessageContext(self)

    def __handle_input(
        self,
        input: InteractiveInput
    ):
        try:
            message = self.__spec.parse_message(
                ParseMessageArgs()
            )
            self.__spec.on_message(
                OnMessageArgs(
                    context=self.__on_message_context,
                    message=message
                )
            )
        except Exception as exn:
            self.__write_error(
                ErrorCodes.UnknownError,
                exn,
                uids = [input.uid]
            )

    def __write_error(
        self,
        error_code: ErrorCodes,
        exn: Exception,
        payload: Optional[str] = None,
        uids: Optional[List[int]] = None 
    ):
        raise NotImplementedError()

    def run(self):
        while(True):
            payload = self.__input_stream.readline()
            try:
                json_value = json.loads(payload)
                value = InteractiveInput(**json_value)
                self.__handle_input(value)
            except json.JSONDecodeError as exn:
                self.__write_error(
                    ErrorCodes.JsonDecodeError,
                    exn,
                    payload,
                )
            except ValidationError as exn:
                self.__write_error(
                    ErrorCodes.PydanticValidationError,
                    exn,
                    payload,
                )
            