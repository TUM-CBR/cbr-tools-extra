import asyncio
from Bio import Entrez as entrez
from Bio.Entrez.Parser import DictionaryElement, ListElement
from Bio.Blast.Record import Blast
from Bio.Blast import NCBIWWW as blast
from Bio.Blast import NCBIXML as blastxml
from concurrent.futures import Executor
from itertools import groupby
import re
from typing import Dict, Iterable, List

from .data import *


K_ID_LIST = 'IDList'
TAXID_REGEX = re.compile(r".+\(taxid:\d+\)\s*")
K_SCIENTIFIC_NAME = "ScientificName"
K_TAX_ID = "TaxId"


def get_tax_ids(organisms : Iterable[str]) -> Iterable[str]:

    organisms = list(organisms)
    query = "+OR+".join(f"\"{organism}\"[All Names]" for organism in organisms)

    search_result = entrez.read(
        entrez.esearch(db="taxonomy", term=query, retmax=len(organisms))
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
        entrez.esearch(db="nucleotide", term=query, retmax=len(accessions))
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


async def find_organisms(
    args : FindOrganismsArgs,
    executor : Executor
) -> CascadeStepResult:
    included_organisms = args.included_organisms
    excluded_organisms = args.excluded_organisms
    num_results = args.num_results
    step = args.step

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

    def query_step(fasta: str) -> Blast:
        buffer = \
            blast.qblast(
                database = "nr",
                sequence = fasta,
                program = "tblastn",
                hitlist_size = num_results,
                alignments = num_results,
                descriptions = num_results,
                max_num_seq = num_results,
                organisms=organisms
            )

        result = blastxml.read(buffer)
        assert isinstance(result, Blast)
        return result

    blast_results = await asyncio.gather(*
        [
            asyncio.get_event_loop().run_in_executor(executor, query_step, fasta)
            for fasta in step.fasta
        ]
    )

    accession_to_organism = get_organisms_for_accesions(
        alignment.accession
        for result in blast_results
        for alignment in result.alignments
    )

    organisms = {}

    for result_ix, result in enumerate(blast_results):
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

    return CascadeStepResult(
        step = step,
        organisms=list(organisms.values())
    )

