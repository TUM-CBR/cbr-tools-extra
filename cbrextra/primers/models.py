from sqlalchemy import create_engine, Column, Float, Integer, String, Sequence
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class PrimerModel(Base):
    __tablename__ = 'primers'
    id = Column(Integer, Sequence('primier_id'), primary_key=True)
    primer_left = Column(String(50))
    primer_right = Column(String(50))
    tm_left = Column(Float)
    tm_right = Column(Float)
    innter_seq = Column(String(50))
    tm_all = Column(Float)
