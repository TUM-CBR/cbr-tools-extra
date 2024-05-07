from Bio import SeqIO
from os import path
from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session, sessionmaker
from typing import Any, Optional, TextIO, cast, List, NamedTuple, Sequence

from .data import DnaSeq, SequenceLoadException
from .models import Base, DnaSeqModel, LogEntryModel, SeqDataFolderModel

class SessionInstance:

    def __init__(self, session: Session) -> None:
        self.__session = session

    def __enter__(self) -> 'SessionInstance':
        self.__session.__enter__()
        return self

    def __exit__(self, *args: Any, **kwargs: Any):
        self.__session.commit()
        self.__session.__exit__(*args, **kwargs)

    def get_paths(self) -> Sequence[str]:

        return cast(
            Sequence[str],
            [
                data_folder.location
                for (data_folder,) in self.__session.execute(select(SeqDataFolderModel))
            ]
        )

    def update_search_paths(self, paths: Sequence[str]) -> None:

        existing = self.get_paths()

        for path in paths:

            if path in existing:
                continue

            folder = SeqDataFolderModel(location = path)
            self.__session.add(folder)

    def update_sequences(self, *seqs: DnaSeq) -> Sequence[DnaSeqModel]:

        result: List[DnaSeqModel] = []

        for seq in seqs:

            query = select(DnaSeqModel).where(DnaSeqModel.seq_id == seq.seq.id)
            existing = [item for (item,) in self.__session.execute(query)]

            if len(existing) > 0:
                model = existing[0]
            else:
                model = DnaSeqModel.from_dna(seq)
                self.__session.add(model)

            result.append(model)

        return result
    
    def add_error(self, error: SequenceLoadException) -> LogEntryModel:
        exn = LogEntryModel.from_exception(error)
        self.__session.add(exn)
        return exn
        
    def as_fasta(self, out_stream: TextIO):

        seqs = (
            model.as_seq_record()
            for (model,) in self.__session.execute(select(DnaSeqModel))
            for model in [cast(DnaSeqModel, model)]
        )

        SeqIO.write(seqs, out_stream, format='fasta')

    def dnaseq_by_id(self, uid: int) -> Optional[DnaSeqModel]:
        query = select(DnaSeqModel).where(DnaSeqModel.id == uid)
        return next(
            (record for (record,) in self.__session.execute(query)),
            None
        )

class SessionManager(NamedTuple):

    session_factory: sessionmaker[Any]

    def session(self) -> SessionInstance:
        return SessionInstance(self.session_factory())

    @classmethod
    def from_file(cls, db_file: str) -> 'SessionManager':

        engine = create_engine(f"sqlite:///{db_file}")

        if not path.exists(db_file):
            Base.metadata.create_all(engine)

        return SessionManager(
            session_factory=sessionmaker(bind=engine)
        )

