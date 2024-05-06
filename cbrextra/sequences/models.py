from sqlalchemy import Column, Integer, Sequence, String
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class DnaSeq(Base):
    __tablename__ = "dnaseq"
    id = Column(Integer, Sequence("dnaseq_id"), primary_key=True)
    tax_id = Column(Integer, nullable=True)
    seq_id = Column(String)
    seq = Column(String)
    seq_file = Column(String)

class LogEntry(Base):
    __tablename__ = "logentry"
    id = Column(Integer, Sequence("logentry_id"), primary_key=True)
    message = Column(String)
    msg_file = Column(String)

class SeqDataFolder(Base):
    __tablename__ = "seqdatafolder"
    id = Column(Integer, Sequence("seqdatafolder_id"), primary_key=True)
    location = Column(String)