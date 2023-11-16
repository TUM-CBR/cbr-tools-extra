from sqlalchemy import Column, Float, ForeignKey, Integer, String, Sequence, Text, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from typing import Any, Dict, Iterable

Base = declarative_base()

class OrganismModel(Base):
    __tablename__ = "organisms"
    tax_id = Column(Integer, primary_key=True)
    name = Column(String(128))
    cascade_organisms = relationship('CascadeOrganism', back_populates='organism')

class CascadeStepOrganismModel(Base):
    __tablename__ = "cascade_step_organism"
    id = Column(Integer, Sequence("cascade_organism_id"), primary_key=True)
    identity = Column(Float)
    sequence_match = String(2**13)
    organism_id = Column(Integer, ForeignKey("organisms.tax_id"))
    organism = relationship("OrganismModel", back_populates='cascade_organisms')
    step = Column(Integer, ForeignKey='cascade_steps.id')