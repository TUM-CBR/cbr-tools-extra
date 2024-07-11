from Bio import Entrez
from Bio.Entrez.Parser import DictionaryElement, ListElement 
from typing import Iterable, List, Optional, Union

from .data import SearchArg, SearchArgs, SearchResult

K_SCIENTIFIC_NAME = "ScientificName"
K_TAX_DB = "taxonomy"
K_GENOME_DB = ""

class SearchResultBuilder:

    def __init__(self, search_id):
        self.__search_id = search_id
        self.__errors = []
    def add_errors(
        self,
        *values: Union[str, Exception]
    ):
        self.__errors += [
            str(value)
            for value in values
        ]

    def build(self) -> SearchResult:
        return SearchResult(
            search_id = self.__search_id,
            errors = self.__errors
        )
class Search:

    def get_names(
        self,
        builder: SearchResultBuilder,
        tax_ids: List[int]
    ) -> List[str]:

        if len(tax_ids) == 0:
            return []

        result = Entrez.read(Entrez.esummary(db=K_TAX_DB, id=tax_ids))
    
        if not isinstance(result, ListElement):
            builder.add_errors(f"Querying taxids failed: {str(result)}")
            return []

        return [
            org[K_SCIENTIFIC_NAME]
            for org in result
        ]

    def run_single(
        self,
        arg: SearchArg
    ):

        builder = SearchResultBuilder(arg.search_id)
        organisms = list(arg.names or [])
        if arg.tax_ids is not None:
            organisms += self.get_names(builder, arg.tax_ids) 

        query = \
            " OR ".join(f"\"{organism}\"[Organism]" for organism in organisms)

        if arg.accession is not None:
            query += " OR " + " OR ".join(f"\"{arg.accession}\"[{field}]" for field in ["PRJA", "AACC", "ACCN"])
        
        result_ids = Entrez.read(Entrez.esearch(db=K_GENOME_DB, term=query))
        
        if not isinstance(result_ids, DictionaryElement):
            builder.add_errors(f"Got unexpected result while quering the ids: {str(result_ids)}")

    def run(
        self,
        args: SearchArgs
    ) -> Iterable[SearchResult]:
        return args
