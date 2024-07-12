# coding: utf-8

"""
    NCBI Datasets API

    ### NCBI Datasets is a resource that lets you easily gather data from NCBI. The Datasets version 2 API is still in alpha, and we're updating it often to add new functionality, iron out bugs and enhance usability. For some larger downloads, you may want to download a [dehydrated zip archive](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download/), and retrieve the individual data files at a later time. 

    The version of the OpenAPI document: v2alpha
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


from __future__ import annotations
import pprint
import re  # noqa: F401
import json

from pydantic import BaseModel, ConfigDict, StrictBool, StrictInt, StrictStr
from typing import Any, ClassVar, Dict, List, Optional
from cbrextra.ncbi.openapi.models.v2_include_tabular_header import V2IncludeTabularHeader
from typing import Optional, Set
from typing_extensions import Self

class V2AssemblySequenceReportsRequest(BaseModel):
    """
    V2AssemblySequenceReportsRequest
    """ # noqa: E501
    accession: Optional[StrictStr] = None
    chromosomes: Optional[List[StrictStr]] = None
    role_filters: Optional[List[StrictStr]] = None
    table_fields: Optional[List[StrictStr]] = None
    count_assembly_unplaced: Optional[StrictBool] = None
    page_size: Optional[StrictInt] = None
    page_token: Optional[StrictStr] = None
    include_tabular_header: Optional[V2IncludeTabularHeader] = None
    table_format: Optional[StrictStr] = None
    __properties: ClassVar[List[str]] = ["accession", "chromosomes", "role_filters", "table_fields", "count_assembly_unplaced", "page_size", "page_token", "include_tabular_header", "table_format"]

    model_config = ConfigDict(
        populate_by_name=True,
        validate_assignment=True,
        protected_namespaces=(),
    )


    def to_str(self) -> str:
        """Returns the string representation of the model using alias"""
        return pprint.pformat(self.model_dump(by_alias=True))

    def to_json(self) -> str:
        """Returns the JSON representation of the model using alias"""
        # TODO: pydantic v2: use .model_dump_json(by_alias=True, exclude_unset=True) instead
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, json_str: str) -> Optional[Self]:
        """Create an instance of V2AssemblySequenceReportsRequest from a JSON string"""
        return cls.from_dict(json.loads(json_str))

    def to_dict(self) -> Dict[str, Any]:
        """Return the dictionary representation of the model using alias.

        This has the following differences from calling pydantic's
        `self.model_dump(by_alias=True)`:

        * `None` is only added to the output dict for nullable fields that
          were set at model initialization. Other fields with value `None`
          are ignored.
        """
        excluded_fields: Set[str] = set([
        ])

        _dict = self.model_dump(
            by_alias=True,
            exclude=excluded_fields,
            exclude_none=True,
        )
        return _dict

    @classmethod
    def from_dict(cls, obj: Optional[Dict[str, Any]]) -> Optional[Self]:
        """Create an instance of V2AssemblySequenceReportsRequest from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "accession": obj.get("accession"),
            "chromosomes": obj.get("chromosomes"),
            "role_filters": obj.get("role_filters"),
            "table_fields": obj.get("table_fields"),
            "count_assembly_unplaced": obj.get("count_assembly_unplaced"),
            "page_size": obj.get("page_size"),
            "page_token": obj.get("page_token"),
            "include_tabular_header": obj.get("include_tabular_header"),
            "table_format": obj.get("table_format")
        })
        return _obj

