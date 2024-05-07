from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sqlalchemy import Column, Integer, Sequence, String
from sqlalchemy.ext.declarative import declarative_base
from typing import cast

from .data import DnaSeq, SequenceLoadException

Base = declarative_base()

class DnaSeqModel(Base):
    __tablename__ = "dnaseq"
    id = Column(Integer, Sequence("dnaseq_id"), primary_key=True)
    tax_id = Column(Integer, nullable=True)
    seq_id = Column(String)
    seq = Column(String)
    seq_file = Column(String)

    @classmethod
    def from_dna(cls, dna: DnaSeq) -> 'DnaSeqModel':
        return DnaSeqModel(
            tax_id = dna.tax_id,
            seq_id = dna.seq.id,
            seq = dna.seq._seq._data,
            seq_file = dna.seq_file
        )
    
    def as_seq_record(self):
        return SeqRecord(
            seq = Seq(cast(str, self.seq)),
            id = str(self.id)
        )

class LogEntryModel(Base):
    __tablename__ = "logentry"
    id = Column(Integer, Sequence("logentry_id"), primary_key=True)
    message = Column(String)
    offending_file = Column(String)

    @classmethod
    def from_exception(cls, exn: SequenceLoadException) -> 'LogEntryModel':
        return LogEntryModel(
            message = exn.message,
            offending_file = exn.file
        )

class SeqDataFolderModel(Base):
    __tablename__ = "seqdatafolder"
    id = Column(Integer, Sequence("seqdatafolder_id"), primary_key=True)
    location = Column(String)