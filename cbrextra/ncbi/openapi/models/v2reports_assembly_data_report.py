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
from cbrextra.ncbi.openapi.models.v2reports_annotation_info import V2reportsAnnotationInfo
from cbrextra.ncbi.openapi.models.v2reports_assembly_info import V2reportsAssemblyInfo
from cbrextra.ncbi.openapi.models.v2reports_assembly_stats import V2reportsAssemblyStats
from cbrextra.ncbi.openapi.models.v2reports_average_nucleotide_identity import V2reportsAverageNucleotideIdentity
from cbrextra.ncbi.openapi.models.v2reports_check_m import V2reportsCheckM
from cbrextra.ncbi.openapi.models.v2reports_organelle_info import V2reportsOrganelleInfo
from cbrextra.ncbi.openapi.models.v2reports_organism import V2reportsOrganism
from cbrextra.ncbi.openapi.models.v2reports_source_database import V2reportsSourceDatabase
from cbrextra.ncbi.openapi.models.v2reports_type_material import V2reportsTypeMaterial
from cbrextra.ncbi.openapi.models.v2reports_wgs_info import V2reportsWGSInfo
from typing import Optional, Set
from typing_extensions import Self

class V2reportsAssemblyDataReport(BaseModel):
    """
    V2reportsAssemblyDataReport
    """ # noqa: E501
    accession: Optional[StrictStr] = None
    current_accession: Optional[StrictStr] = None
    paired_accession: Optional[StrictStr] = None
    source_database: Optional[V2reportsSourceDatabase] = None
    organism: Optional[V2reportsOrganism] = None
    assembly_info: Optional[V2reportsAssemblyInfo] = None
    assembly_stats: Optional[V2reportsAssemblyStats] = None
    organelle_info: Optional[List[V2reportsOrganelleInfo]] = None
    annotation_info: Optional[V2reportsAnnotationInfo] = None
    wgs_info: Optional[V2reportsWGSInfo] = None
    type_material: Optional[V2reportsTypeMaterial] = None
    checkm_info: Optional[V2reportsCheckM] = None
    average_nucleotide_identity: Optional[V2reportsAverageNucleotideIdentity] = None
    __properties: ClassVar[List[str]] = ["accession", "current_accession", "paired_accession", "source_database", "organism", "assembly_info", "assembly_stats", "organelle_info", "annotation_info", "wgs_info", "type_material", "checkm_info", "average_nucleotide_identity"]

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
        """Create an instance of V2reportsAssemblyDataReport from a JSON string"""
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
        # override the default output from pydantic by calling `to_dict()` of organism
        if self.organism:
            _dict['organism'] = self.organism.to_dict()
        # override the default output from pydantic by calling `to_dict()` of assembly_info
        if self.assembly_info:
            _dict['assembly_info'] = self.assembly_info.to_dict()
        # override the default output from pydantic by calling `to_dict()` of assembly_stats
        if self.assembly_stats:
            _dict['assembly_stats'] = self.assembly_stats.to_dict()
        # override the default output from pydantic by calling `to_dict()` of each item in organelle_info (list)
        _items = []
        if self.organelle_info:
            for _item in self.organelle_info:
                if _item:
                    _items.append(_item.to_dict())
            _dict['organelle_info'] = _items
        # override the default output from pydantic by calling `to_dict()` of annotation_info
        if self.annotation_info:
            _dict['annotation_info'] = self.annotation_info.to_dict()
        # override the default output from pydantic by calling `to_dict()` of wgs_info
        if self.wgs_info:
            _dict['wgs_info'] = self.wgs_info.to_dict()
        # override the default output from pydantic by calling `to_dict()` of type_material
        if self.type_material:
            _dict['type_material'] = self.type_material.to_dict()
        # override the default output from pydantic by calling `to_dict()` of checkm_info
        if self.checkm_info:
            _dict['checkm_info'] = self.checkm_info.to_dict()
        # override the default output from pydantic by calling `to_dict()` of average_nucleotide_identity
        if self.average_nucleotide_identity:
            _dict['average_nucleotide_identity'] = self.average_nucleotide_identity.to_dict()
        return _dict

    @classmethod
    def from_dict(cls, obj: Optional[Dict[str, Any]]) -> Optional[Self]:
        """Create an instance of V2reportsAssemblyDataReport from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "accession": obj.get("accession"),
            "current_accession": obj.get("current_accession"),
            "paired_accession": obj.get("paired_accession"),
            "source_database": obj.get("source_database"),
            "organism": V2reportsOrganism.from_dict(obj["organism"]) if obj.get("organism") is not None else None,
            "assembly_info": V2reportsAssemblyInfo.from_dict(obj["assembly_info"]) if obj.get("assembly_info") is not None else None,
            "assembly_stats": V2reportsAssemblyStats.from_dict(obj["assembly_stats"]) if obj.get("assembly_stats") is not None else None,
            "organelle_info": [V2reportsOrganelleInfo.from_dict(_item) for _item in obj["organelle_info"]] if obj.get("organelle_info") is not None else None,
            "annotation_info": V2reportsAnnotationInfo.from_dict(obj["annotation_info"]) if obj.get("annotation_info") is not None else None,
            "wgs_info": V2reportsWGSInfo.from_dict(obj["wgs_info"]) if obj.get("wgs_info") is not None else None,
            "type_material": V2reportsTypeMaterial.from_dict(obj["type_material"]) if obj.get("type_material") is not None else None,
            "checkm_info": V2reportsCheckM.from_dict(obj["checkm_info"]) if obj.get("checkm_info") is not None else None,
            "average_nucleotide_identity": V2reportsAverageNucleotideIdentity.from_dict(obj["average_nucleotide_identity"]) if obj.get("average_nucleotide_identity") is not None else None
        })
        return _obj

