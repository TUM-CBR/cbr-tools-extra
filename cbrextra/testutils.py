import json
from os import path
from pydantic import BaseModel
from typing import Any, IO, List, Type, TypeVar

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