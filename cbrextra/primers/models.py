from sqlalchemy import Column, Float, Integer, String, Sequence
from sqlalchemy.ext.declarative import declarative_base

from .data import PrimerResult

Base = declarative_base()

class PrimerModel(Base):
    __tablename__ = 'primers'
    id = Column(Integer, Sequence('primier_id'), primary_key=True)
    primer_left = Column(String(50))
    primer_right = Column(String(50))
    tm_left = Column(Float)
    tm_right = Column(Float)
    inner_seq = Column(String(50))
    tm_all = Column(Float)
    group = Column(Integer)
    codon = Column(String(3))

    @staticmethod
    def from_result(group: int, codon: str, result: PrimerResult) -> 'PrimerModel':
        return PrimerModel(
            primer_left = result.left_primer,
            primer_right = result.right_primer,
            tm_left = result.tm_left,
            tm_right = result.tm_right,
            inner_seq = result.inner_seq,
            tm_all = result.tm_all,
            group = group,
            codon = codon
        )