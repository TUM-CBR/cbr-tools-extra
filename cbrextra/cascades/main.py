from Bio import Entrez as entrez
from Bio.Entrez.Parser import DictionaryElement, ListElement
from Bio.Blast import NCBIWWW as blast
from Bio.Blast import NCBIXML as blastxml
from concurrent.futures import ThreadPoolExecutor
from itertools import groupby
import re
from typing import Any, Dict, Generator, Iterable, Iterator, List, NamedTuple, Optional

from .data import *


K_ID_LIST = 'IDList'
TAXID_REGEX = re.compile(r".+\(taxid:\d+\)\s*")
K_SCIENTIFIC_NAME = "ScientificName"
K_TAX_ID = "TaxId"


def get_tax_ids(organisms : Iterable[str]) -> Iterable[str]:

    query = "+OR+".join(f"\"{organism}\"[All Names]" for organism in organisms)

    search_result = entrez.read(
        entrez.esearch(db="taxonomy", term=query)
    )

    if search_result is None:
        raise Exception(f"Query failed: {query}")

    assert isinstance(search_result, DictionaryElement)

    ids = search_result[K_ID_LIST]
    summary = entrez.read(entrez.esummary(db="taxonomy", id=ids))

    if summary is None:
        raise Exception(f"Failed to query ids: {', '.join(ids)}")

    assert isinstance(summary, ListElement)

    return [
        f"{item[K_SCIENTIFIC_NAME]} (taxid:{item[K_TAX_ID]})"
        for item in summary
    ]


def ensure_taxid(organisms: List[str]) -> List[str]:

    groups = groupby(organisms, lambda o: True if TAXID_REGEX.match(o) else False)

    return [
        organism
        for has_taxid, classified in groups
        for organism in (classified if has_taxid else get_tax_ids(classified))
    ]


K_ID_LIST = "IdList"
K_ACCESSION_VERSION = "AccessionVersion"
K_SCIENTIFIC_NAME = "ScientificName"


def get_organisms_for_accesions(accessions: Iterable[str]) -> Dict[str, Organism]:
    """Given a list of accession names, it will search in the entrez nucleotide database
    for the taxonomy id of the organism of that accession. It will then use the id to
    search the taxonomy database to recover the organism's scientific name"""

    accessions = set(accessions)
    query = " OR ".join(
        f"\"{accession}\"[Accession]"
        for accession in accessions
    )
    articles_iter = entrez.read(
        entrez.esearch(db="nucleotide", term=query)
    )

    ids = []

    for articles in [articles_iter]:
        assert(isinstance(articles, DictionaryElement))
        ids += articles[K_ID_LIST]

    article_summaries = entrez.parse(
        entrez.esummary(
            db = "nucleotide",
            id=ids
        )
    )

    accession_to_tid = {}
    for summary in article_summaries:
        assert(isinstance(summary, DictionaryElement))
        accession_with_version = summary[K_ACCESSION_VERSION].split(".")
        accession = accession_with_version[0]
        version = accession_with_version[1] if len(accession_with_version) == 2 else None
        tid = int(summary[K_TAX_ID])

        accession_to_tid[accession] = tid

        if version is not None:
            accession_to_tid[f"{accession}.{version}"] = tid

    taxids = list(set(accession_to_tid.values()))
    taxonomy_summary = entrez.parse(
        entrez.esummary(
            db = "taxonomy",
            id = taxids
        )
    )

    taxid_to_organism_entry = {}
    for entry in taxonomy_summary:
        assert(isinstance(entry, DictionaryElement))
        taxid = int(entry[K_TAX_ID])
        taxid_to_organism_entry[taxid] = entry

    return dict(
        (accession, Organism(tax_id=tax_id, name = organism_entry[K_SCIENTIFIC_NAME]))
        for accession,tax_id in accession_to_tid.items()
        for organism_entry in [taxid_to_organism_entry[tax_id]]
    )


default_include = [
    "Bacteria (taxid:2)",
    #"Archaea (taxid:2157)"
]
default_exclude = ["synthetic constructs (taxid:32630)"]


def find_cascades(
    cascade : List[CascadeStep],
    excluded_organisms: Optional[List[str]] = None,
    included_organisms: Optional[List[str]] = None,
    num_results = 5000
) -> CascadeReesult:

    if included_organisms is None:
        included_organisms = list(default_include)

    if excluded_organisms is None:
        excluded_organisms = list(default_exclude)

    included_organisms = ensure_taxid(included_organisms)
    excluded_organisms = ensure_taxid(excluded_organisms)

    organisms = {}

    for organism in included_organisms:
        organisms[organism] = False

    for organism in excluded_organisms:
        organisms[organism] = True

    def query_step(step : CascadeStep) -> List[Any]:
        buffers = [
            blast.qblast(
                database = "nr",
                sequence = fasta,
                program = "tblastn",
                hitlist_size = num_results,
                alignments = num_results,
                descriptions = num_results,
                organisms=organisms
            )
            for fasta in step.fasta
        ]

        return [
            blastxml.read(buffer)
            for buffer in buffers
        ]

    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = [
            executor.submit(
                query_step,
                step
            )
            for step in cascade
        ]

        #Todo: maybe not keep everything in memory
        blast_results = [
            future.result()
            for future in futures
        ]

    accession_to_organism = get_organisms_for_accesions(
        alignment.accession
        for results in blast_results
        for result in results
        for alignment in result.alignments
    )

    cascade_steps = []

    for i, results in enumerate(blast_results):

        step = cascade[i]
        organisms = {}

        for result_ix, result in results:
            sequence = step.sequences[result_ix]
            for alignment in result.alignments:
                accession = alignment.accession
                organism = accession_to_organism[accession]
                closets_match = max(
                    alignment.hsps,
                    key = lambda hsps: hsps.identities
                )
                cascade_organism = CascadeStepOrganism(
                    organism=organism,
                    identity=closets_match.identities/len(sequence),
                    sequence_match=closets_match.sbjct
                )
                update_if_better(organisms, cascade_organism)

        cascade_steps.append(
            CascadeStepResult(
                step = step,
                organisms = list(organisms.values())
            )
        )

    return CascadeReesult(
        steps = cascade_steps
    )

