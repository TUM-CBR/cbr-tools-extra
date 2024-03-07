from io import StringIO
import json
from os import path
from pydantic import BaseModel
from typing import Any, Dict, IO, Generic, Iterable, List, Optional, Type, TypeVar

from .core.interactive import InteractiveSpec, run_interactive
from .core.interactive.data import EntityType, InteractiveError, InteractiveInput, InteractiveOutput

def get_resource_path(
    file_loc: str,
    *resoruce: str
) -> str:
    return path.join(
        path.dirname(file_loc),
        "resources",
        *resoruce
    )

def open_resource(
    file_loc: str,
    *resource: str,
    mode: str = 'r'
) -> IO[Any]:
    return open(get_resource_path(file_loc, *resource), mode)

TModel = TypeVar('TModel', bound=BaseModel)

def open_models(
    _class: Type[TModel],
    file_loc: str,
    *resoruce: str
) -> List[TModel]:
    
    with open_resource(file_loc, *resoruce) as json_model:
        return [
            _class(**model)
            for model in json.load(json_model)
        ]
    
TMessageModelIn = TypeVar('TMessageModelIn', bound=BaseModel)
TMessageModelOut = TypeVar('TMessageModelOut', bound=BaseModel)

class TestInteractiveOutput(Generic[TMessageModelIn, TMessageModelOut]):
    in_messages: List[TMessageModelIn]
    result: Optional[TMessageModelOut]
    error: Optional[InteractiveError]

    @classmethod
    def from_output(
        cls,
        out_json: Dict[Any, Any],
        in_messages: Dict[int, TMessageModelIn],
        out_model: Type[TMessageModelOut]
    ) -> 'TestInteractiveOutput[TMessageModelIn, TMessageModelOut]':

        result: TestInteractiveOutput[TMessageModelIn, TMessageModelOut] = TestInteractiveOutput()
        i_output = InteractiveOutput(**out_json)
        result.error = i_output.error

        value = i_output.value
        if value is not None:
            result.in_messages = [
                in_messages[mid]
                for mid in value.input_uids
            ]
            result.result = out_model(**value.payload)
        else:
            result.in_messages = []
            result.result = None

        return result

def test_interactive(
    spec: InteractiveSpec[TMessageModelIn, TMessageModelOut],
    input: Iterable[TMessageModelIn],
    _in_model: Type[TMessageModelIn],
    out_model: Type[TMessageModelOut]
) -> Iterable[TestInteractiveOutput[TMessageModelIn, TMessageModelOut]]:
    messages: Dict[int, TMessageModelIn] = {}
    
    with StringIO() as input_stream, \
        StringIO() as output_stream:

        for uid, message in enumerate(input):
            messages[uid] = message
            model_json = message.model_dump()
            in_message = InteractiveInput(uid=uid, entity_type=EntityType.MESSAGE.value, payload=model_json)
            json.dump(in_message.model_dump(), input_stream)
            input_stream.write("\n")

        input_stream.seek(0)

        run_interactive(spec, input_stream=input_stream, output_stream=output_stream)

        output_stream.seek(0)

        for line in output_stream.readlines():
            out_json = json.loads(line)
            yield TestInteractiveOutput.from_output(
                out_json,
                messages,
                out_model
            )