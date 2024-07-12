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
from cbrextra.ncbi.openapi.models.v2reports_ani_match import V2reportsANIMatch
from cbrextra.ncbi.openapi.models.v2reports_ani_type_category import V2reportsANITypeCategory
from cbrextra.ncbi.openapi.models.v2reports_average_nucleotide_identity_match_status import V2reportsAverageNucleotideIdentityMatchStatus
from cbrextra.ncbi.openapi.models.v2reports_average_nucleotide_identity_taxonomy_check_status import V2reportsAverageNucleotideIdentityTaxonomyCheckStatus
from typing import Optional, Set
from typing_extensions import Self

class V2reportsAverageNucleotideIdentity(BaseModel):
    """
    V2reportsAverageNucleotideIdentity
    """ # noqa: E501
    taxonomy_check_status: Optional[V2reportsAverageNucleotideIdentityTaxonomyCheckStatus] = None
    match_status: Optional[V2reportsAverageNucleotideIdentityMatchStatus] = None
    submitted_organism: Optional[StrictStr] = None
    submitted_species: Optional[StrictStr] = None
    category: Optional[V2reportsANITypeCategory] = None
    submitted_ani_match: Optional[V2reportsANIMatch] = None
    best_ani_match: Optional[V2reportsANIMatch] = None
    comment: Optional[StrictStr] = None
    __properties: ClassVar[List[str]] = ["taxonomy_check_status", "match_status", "submitted_organism", "submitted_species", "category", "submitted_ani_match", "best_ani_match", "comment"]

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
        """Create an instance of V2reportsAverageNucleotideIdentity from a JSON string"""
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
        # override the default output from pydantic by calling `to_dict()` of submitted_ani_match
        if self.submitted_ani_match:
            _dict['submitted_ani_match'] = self.submitted_ani_match.to_dict()
        # override the default output from pydantic by calling `to_dict()` of best_ani_match
        if self.best_ani_match:
            _dict['best_ani_match'] = self.best_ani_match.to_dict()
        return _dict

    @classmethod
    def from_dict(cls, obj: Optional[Dict[str, Any]]) -> Optional[Self]:
        """Create an instance of V2reportsAverageNucleotideIdentity from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "taxonomy_check_status": obj.get("taxonomy_check_status"),
            "match_status": obj.get("match_status"),
            "submitted_organism": obj.get("submitted_organism"),
            "submitted_species": obj.get("submitted_species"),
            "category": obj.get("category"),
            "submitted_ani_match": V2reportsANIMatch.from_dict(obj["submitted_ani_match"]) if obj.get("submitted_ani_match") is not None else None,
            "best_ani_match": V2reportsANIMatch.from_dict(obj["best_ani_match"]) if obj.get("best_ani_match") is not None else None,
            "comment": obj.get("comment")
        })
        return _obj

