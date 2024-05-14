from io import StringIO
import json
from os import path
import tempfile
from typing import Any, Optional, cast, Dict, Iterable, List, Sequence, TextIO
import warnings

from .blast import Blast
from .clonemanager import CMLoaderFactory
from .data import DnaSeq, QueryResult, QueryResults, SequencesContext, SequenceLoadException
from .excelloader import ExcelLoaderFactory
from .manager import SessionInstance, SessionManager
from .sequence import SeqLoaderItem, SeqLoaderManager

def load_sequences(paths: Sequence[str]) -> Iterable[SeqLoaderItem]:
    loader = SeqLoaderManager(
        paths,
        [
            CMLoaderFactory(),
            ExcelLoaderFactory()
        ]
    )

    return loader.load_files()

def run_scan(context: SequencesContext, extra_folders: Sequence[str]):
    db_file = context.db_file

    manager = SessionManager.from_file(db_file)

    with manager.session() as session:
        session.update_search_paths(extra_folders)
        extra_folders = session.get_paths()

    # First we update the sequences database by scanning all
    # the sequence files (in various formats)
    for sequence in load_sequences(extra_folders):

        with manager.session() as session:
            if isinstance(sequence, DnaSeq):
                session.update_sequences(sequence)
            else:
                assert isinstance(sequence, SequenceLoadException)
                session.add_error(sequence)

    with tempfile.TemporaryDirectory() as workdir:

        blast = Blast(context.blast_context)
        fasta_db = path.join(workdir, "all.fasta")

        with open(fasta_db, 'w') as fasta_stream \
            , manager.session() as session:
            session.as_fasta(fasta_stream)

        blast.update_dbs(
            context.db_blast_name,
            fasta_db
        )

K_BLAST_OUTPUT = "BlastOutput2"
K_REPORT = "report"
K_RESULTS = "results"
K_SEARCH = "search"
K_QUERY_ID = "query_id"
K_QUERY_LEN = "query_len"
K_HITS = "hits"
K_DESCRIPTION = "description"
K_ACCESSION = "accession"
K_HSPS = "hsps"
K_IDENTITY = "identity"
K_QUERY_SEQ = "qseq"
K_QUERY_HSEQ = "hseq"
K_MIDLINE = "midline"

def parse_blast_results(
    session: SessionInstance,
    result: Dict[Any, Any]
) -> Sequence[QueryResult]:
    
    output_region: List[Dict[Any, Any]] = result[K_BLAST_OUTPUT]
    return_value: List[QueryResult] = []

    for output in output_region:
        report: Dict[Any, Any] = output.get(K_REPORT, {})
        results: Dict[Any, Any] = report.get(K_RESULTS, {})
        search: Dict[Any, Any] = results.get(K_SEARCH, {})
        query_id: Optional[str] = search.get(K_QUERY_ID)
        query_len: Optional[int] = search.get(K_QUERY_LEN)
        hits: Optional[List[Dict[Any, Any]]] = search.get(K_HITS)

        if query_id is None or query_len is None or hits is None:
            raise Exception(
                f"Unable to parse blast response: {output}"
            )

        for hit in hits:
            descriptions: Optional[List[Dict[Any, Any]]] = hit.get(K_DESCRIPTION)
            hsps: Optional[List[Dict[Any, Any]]] = hit[K_HSPS]

            if descriptions is None or hsps is None:
                warnings.warn(
                    f"Unable to parse hit {hit}"
                )
                continue

            for (description, hsp) in zip(descriptions, hsps):

                accession_str: Optional[str] = description.get(K_ACCESSION)
                identity_count: Optional[int] = hsp.get(K_IDENTITY)
                query_sequence: Optional[str] = hsp.get(K_QUERY_SEQ)
                hit_sequence: Optional[str] = hsp.get(K_QUERY_HSEQ)
                mid_line: Optional[str] = hsp.get(K_MIDLINE)

                if accession_str is None or identity_count is None \
                    or query_sequence is None or hit_sequence is None \
                    or mid_line is None:
                    warnings.warn(
                        f"Unable to parse {description}, {hit}"
                    )
                    continue

                accession = int(accession_str)
                model = session.dnaseq_by_id(accession)
                if model is None:
                    warnings.warn(
                        f"The sequence {accession} is not in a file."
                    )
                    continue

                return_value.append(
                    QueryResult(
                        id = str(model.seq_id),
                        query_id=query_id,
                        organism=None,
                        identity=identity_count/query_len,
                        query_sequence=query_sequence,
                        hit_sequence=hit_sequence,
                        mid_line=mid_line,
                        file_location=cast(str, model.seq_file)
                    )
                )

    return return_value


def run_query_tblastn(context: SequencesContext, query: TextIO, result_stream: TextIO):

    blast = Blast(context.blast_context)

    with StringIO() as out_stream:
        blast.query_tblastn(
            context.db_blast_name,
            query,
            out_stream
        )

        out_stream.seek(0)
        result = json.load(out_stream)

    db_file = context.db_file

    manager = SessionManager.from_file(db_file)

    with manager.session() as session:
        results = parse_blast_results(session, result)

    result_stream.write(
        QueryResults(results=results).model_dump_json()
    )