# coding: utf-8

"""
    NCBI Datasets API

    ### NCBI Datasets is a resource that lets you easily gather data from NCBI. The Datasets version 2 API is still in alpha, and we're updating it often to add new functionality, iron out bugs and enhance usability. For some larger downloads, you may want to download a [dehydrated zip archive](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download/), and retrieve the individual data files at a later time. 

    The version of the OpenAPI document: v2alpha
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


from __future__ import annotations
import json
from enum import Enum
from typing_extensions import Self


class V2VirusTableField(str, Enum):
    """
    V2VirusTableField
    """

    """
    allowed enum values
    """
    UNSPECIFIED = 'unspecified'
    NUCLEOTIDE_ACCESSION = 'nucleotide_accession'
    SPECIES_TAX_ID = 'species_tax_id'
    SPECIES_NAME = 'species_name'
    GENUS = 'genus'
    FAMILY = 'family'
    NUCLEOTIDE_LENGTH = 'nucleotide_length'
    ISOLATE_NAME = 'isolate_name'
    SEQUENCE_TYPE = 'sequence_type'
    NUC_COMPLETENESS = 'nuc_completeness'
    GEO_LOCATION = 'geo_location'
    US_STATE = 'us_state'
    HOST_NAME = 'host_name'
    HOST_TAX_ID = 'host_tax_id'
    COLLECTION_DATE = 'collection_date'
    BIOPROJECT = 'bioproject'
    BIOSAMPLE = 'biosample'
    POLYPROTEIN_NAME = 'polyprotein_name'
    PROTEIN_NAME = 'protein_name'
    PROTEIN_ACCESSION = 'protein_accession'
    PROTEIN_SYNONYM = 'protein_synonym'
    CDS_SPAN = 'cds_span'

    @classmethod
    def from_json(cls, json_str: str) -> Self:
        """Create an instance of V2VirusTableField from a JSON string"""
        return cls(json.loads(json_str))

