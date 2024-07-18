from Bio import Entrez
from Bio.Entrez.Parser import DictionaryElement, ListElement 
from typing import Iterable, List, Optional, Union

from ..ncbi.openapi import GenomeApi

from .data import SearchArg, SearchArgs, SearchResult, SearchResultRecord

K_SCIENTIFIC_NAME = "ScientificName"
K_TAX_DB = "taxonomy"
K_GENOME_DB = ""
K_ID_LIST = "IdList"

class SearchResultBuilder:

    def __init__(self, search_id):
        self.__search_id = search_id
        self.__errors = []
        self.__results = {}

    def add_errors(
        self,
        *values: Union[str, Exception]
    ):
        self.__errors += [
            str(value)
            for value in values
        ]

    def add_results(
        self,
        *values: SearchResultRecord
    ):
        for value in values:
            self.__results[value.accession] = value

    def build(self) -> SearchResult:
        return SearchResult(
            search_id = self.__search_id,
            errors = self.__errors,
            records = list(self.__results.values())
        )

class Search:

    def __init__(
        self,
        email: Optional[str] = None,
        api_key: Optional[str] = None
    ):
        self.__genome_api = GenomeApi()

    def get_taxids_from_names(
        self,
        builder: SearchResultBuilder,
        names: List[str]
    ) -> List[int]:

        if len(names) == 0:
            return []

        term = " OR ".join(f'"{name}"[ORGN]' for name in names if len(name) > 0)
        retmax = 10
        result = Entrez.read(Entrez.esearch(db=K_TAX_DB, term = term, retmax = retmax))
    
        if not isinstance(result, DictionaryElement):
            builder.add_errors(f"Querying taxids failed: {str(result)}")
            return []

        taxids = result.get(K_ID_LIST)

        if not isinstance(taxids, ListElement):
            builder.add_errors(f"Querying taxids failed: {str(result)}")
            return []

        return taxids

    def get_records_from_taxids(
        self,
        builder: SearchResultBuilder,
        taxids: List[int]
    ) -> List[SearchResultRecord]:
        results = self.__genome_api.genome_dataset_reports_by_taxon(taxons=[str(taxid) for taxid in taxids]) 

        if results.reports is None:
            builder.add_errors(f"Failed to query genome dataset reports: {results}")
            return []

        return [
            report
            for result in results.reports
            for report in [SearchResultRecord.from_assembly_data_report(result)]
                if report is not None
        ]

    def get_records_from_accessions(
        self,
        builder: SearchResultBuilder,
        accessions: List[str]
    ) -> List[SearchResultRecord]:

        results = self.__genome_api.genome_dataset_report(accessions=accessions)

        if results.reports is None:
            builder.add_errors(f"Failed to query genome dataset reports: {results}")
            return []

        return [
            report
            for result in results.reports
            for report in [SearchResultRecord.from_assembly_data_report(result)]
                if report is not None
        ]

    def run_single(
        self,
        builder: SearchResultBuilder,
        arg: SearchArg
    ): 

        taxids = arg.tax_ids or []

        if arg.names is not None:
            taxids += self.get_taxids_from_names(builder, arg.names)

        records: List[SearchResultRecord] = []
        if len(taxids) > 0:
            records += self.get_records_from_taxids(builder, taxids)

        if arg.accession is not None:
            records += self.get_records_from_accessions(builder, accessions=[arg.accession])

        builder.add_results(*records)

    def run(
        self,
        args: SearchArgs
    ) -> Iterable[SearchResult]:

        for arg in args.searches:
            builder = SearchResultBuilder(arg.search_id)
            try:
                self.run_single(builder, arg)
            except Exception as e:
                builder.add_errors(e)
            finally:
                yield builder.build()


