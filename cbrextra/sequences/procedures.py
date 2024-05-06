from os import path
import tempfile
from typing import Iterable, Sequence

from .blast import Blast
from .clonemanager import CMLoader
from .data import DnaSeq, SequencesContext, SequenceLoadException
from .manager import SessionManager
from .sequence import SeqLoaderItem, SeqLoaderManager

def load_sequences(paths: Sequence[str]) -> Iterable[SeqLoaderItem]:
    loader = SeqLoaderManager(
        paths,
        [CMLoader()]
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

        