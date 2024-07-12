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
from cbrextra.ncbi.openapi.models.v2_image_size import V2ImageSize
from typing import Optional, Set
from typing_extensions import Self

class V2TaxonomyImageMetadataResponse(BaseModel):
    """
    V2TaxonomyImageMetadataResponse
    """ # noqa: E501
    tax_id: Optional[StrictStr] = None
    src: Optional[StrictStr] = None
    license: Optional[StrictStr] = None
    attribution: Optional[StrictStr] = None
    source: Optional[StrictStr] = None
    image_sizes: Optional[List[V2ImageSize]] = None
    format: Optional[StrictStr] = None
    __properties: ClassVar[List[str]] = ["tax_id", "src", "license", "attribution", "source", "image_sizes", "format"]

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
        """Create an instance of V2TaxonomyImageMetadataResponse from a JSON string"""
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
        """Create an instance of V2TaxonomyImageMetadataResponse from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "tax_id": obj.get("tax_id"),
            "src": obj.get("src"),
            "license": obj.get("license"),
            "attribution": obj.get("attribution"),
            "source": obj.get("source"),
            "image_sizes": obj.get("image_sizes"),
            "format": obj.get("format")
        })
        return _obj


