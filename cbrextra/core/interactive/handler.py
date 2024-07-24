import json
from pydantic import ValidationError
import sys
import traceback
from typing import TextIO

from ..atomic import AtomicCounter
from .data import *
from .shared import *

class InteractiveHandler(Generic[TMessageIn, TMessageOut]):

    MESSAGE_ID_COUNTER = AtomicCounter()

    class OnMessageContext(OnMessageContextBase[Any, Any]):

        def __init__(
            self,
            handler: 'InteractiveHandler[TMessageIn, TMessageOut]'
        ) -> None:
            super().__init__()
            self.__handler = handler

        def send(
            self,
            message: TMessageOut,
            message_ids : List[int]
        ) -> None:
            self.__handler.__write_message__(message, message_ids)

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

    def __write_message__(
        self,
        message: TMessageOut,
        message_ids : List[int]
    ):
        json_value = self.__spec.serialize_message(
            SerializeMessageArgs(message)
        )

        output = InteractiveOutput(
            uid = self.MESSAGE_ID_COUNTER.increment(),
            value = InteractiveValue(
                input_uids = message_ids,
                payload = json_value
            ),
            error = None
        )
        self.__write_output(output)

    def __write_output(self, output: InteractiveOutput):
        self.__output_stream.write(
            f"{output.model_dump_json()}\n"
        )
        self.__output_stream.flush()

    def __handle_message_input(
        self,
        input: InteractiveInput,
        payload : Dict[Any, Any]
    ):

        try:
            message = self.__spec.parse_message(
                ParseMessageArgs(
                    payload=payload
                )
            )
            self.__spec.on_message(
                OnMessageArgs(
                    context=self.__on_message_context,
                    message=message,
                    header=InputHeader(
                        uid=input.uid
                    )
                )
            )
        except Exception as exn:
            self.__write_error(
                ErrorCodes.UnknownError,
                traceback.format_exc(),
                uids = [input.uid]
            )

    def __handle_input(
        self,
        input: InteractiveInput
    ) -> bool:

        if input.entity_type == EntityType.STOP.value:
            return False
        elif input.entity_type == EntityType.MESSAGE.value:
            if input.payload is None:
                self.__write_error(
                    ErrorCodes.UnknownInputType,
                    "Payload cannot be None",
                    uids=[input.uid]
                )
            else:
                self.__handle_message_input(
                    input,
                    input.payload
                )
        else:
            self.__write_error(
                ErrorCodes.UnknownInputType,
                "Unknown message input type",
                uids=[input.uid]
            )

        return True

    def __write_error(
        self,
        error_code: ErrorCodes,
        message: str,
        payload: Optional[str] = None,
        uids: Optional[List[int]] = None 
    ):
        output = InteractiveOutput(
            uid = self.MESSAGE_ID_COUNTER.increment(),
            error = InteractiveError(
                input_uids = uids or [],
                error_code = error_code,
                payload = payload,
                message = message
            ),
            value = None
        )

        self.__write_output(output)

    def run(self):
        while(True):
            payload = self.__input_stream.readline()

            if len(payload) == 0:
                return
            try:
                json_value = json.loads(payload)
                value = InteractiveInput(**json_value)
                
                if not self.__handle_input(value):
                    return
            except json.JSONDecodeError as exn:
                self.__write_error(
                    ErrorCodes.JsonDecodeError,
                    traceback.format_exc(),
                    payload,
                )
            except ValidationError as exn:
                self.__write_error(
                    ErrorCodes.PydanticValidationError,
                    traceback.format_exc(),
                    payload,
                )
            
def run_interactive(
    spec: InteractiveSpec[TMessageIn, TMessageOut],
    input_stream : TextIO = sys.stdin,
    output_stream : TextIO = sys.stdout
):
    InteractiveHandler(spec, input_stream, output_stream).run()
