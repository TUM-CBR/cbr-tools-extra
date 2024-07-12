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

from pydantic import BaseModel, ConfigDict, Field, StrictBool, StrictInt, StrictStr
from typing import Any, ClassVar, Dict, List, Optional
from cbrextra.ncbi.openapi.models.v2reports_message import V2reportsMessage
from cbrextra.ncbi.openapi.models.v2reports_organelle import V2reportsOrganelle
from typing import Optional, Set
from typing_extensions import Self

class V2reportsOrganelleDataReports(BaseModel):
    """
    V2reportsOrganelleDataReports
    """ # noqa: E501
    messages: Optional[List[V2reportsMessage]] = None
    reports: Optional[List[V2reportsOrganelle]] = None
    total_count: Optional[StrictInt] = None
    next_page_token: Optional[StrictStr] = None
    report_type: Optional[StrictStr] = Field(default=None, alias="_report_type")
    report_fields: Optional[List[StrictStr]] = Field(default=None, alias="_report_fields")
    first_page: Optional[StrictBool] = Field(default=None, alias="_first_page")
    report_format: Optional[StrictStr] = Field(default=None, alias="_report_format")
    __properties: ClassVar[List[str]] = ["messages", "reports", "total_count", "next_page_token", "_report_type", "_report_fields", "_first_page", "_report_format"]

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
        """Create an instance of V2reportsOrganelleDataReports from a JSON string"""
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
        # override the default output from pydantic by calling `to_dict()` of each item in messages (list)
        _items = []
        if self.messages:
            for _item in self.messages:
                if _item:
                    _items.append(_item.to_dict())
            _dict['messages'] = _items
        # override the default output from pydantic by calling `to_dict()` of each item in reports (list)
        _items = []
        if self.reports:
            for _item in self.reports:
                if _item:
                    _items.append(_item.to_dict())
            _dict['reports'] = _items
        return _dict

    @classmethod
    def from_dict(cls, obj: Optional[Dict[str, Any]]) -> Optional[Self]:
        """Create an instance of V2reportsOrganelleDataReports from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "messages": [V2reportsMessage.from_dict(_item) for _item in obj["messages"]] if obj.get("messages") is not None else None,
            "reports": [V2reportsOrganelle.from_dict(_item) for _item in obj["reports"]] if obj.get("reports") is not None else None,
            "total_count": obj.get("total_count"),
            "next_page_token": obj.get("next_page_token"),
            "_report_type": obj.get("_report_type"),
            "_report_fields": obj.get("_report_fields"),
            "_first_page": obj.get("_first_page"),
            "_report_format": obj.get("_report_format")
        })
        return _obj


