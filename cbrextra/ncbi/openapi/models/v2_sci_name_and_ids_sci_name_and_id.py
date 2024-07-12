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

from pydantic import BaseModel, ConfigDict, StrictStr
from typing import Any, ClassVar, Dict, List, Optional
from cbrextra.ncbi.openapi.models.v2reports_rank_type import V2reportsRankType
from typing import Optional, Set
from typing_extensions import Self

class V2SciNameAndIdsSciNameAndId(BaseModel):
    """
    V2SciNameAndIdsSciNameAndId
    """ # noqa: E501
    sci_name: Optional[StrictStr] = None
    tax_id: Optional[StrictStr] = None
    common_name: Optional[StrictStr] = None
    matched_term: Optional[StrictStr] = None
    rank: Optional[V2reportsRankType] = None
    group_name: Optional[StrictStr] = None
    __properties: ClassVar[List[str]] = ["sci_name", "tax_id", "common_name", "matched_term", "rank", "group_name"]

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
        """Create an instance of V2SciNameAndIdsSciNameAndId from a JSON string"""
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
        """Create an instance of V2SciNameAndIdsSciNameAndId from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "sci_name": obj.get("sci_name"),
            "tax_id": obj.get("tax_id"),
            "common_name": obj.get("common_name"),
            "matched_term": obj.get("matched_term"),
            "rank": obj.get("rank"),
            "group_name": obj.get("group_name")
        })
        return _obj


