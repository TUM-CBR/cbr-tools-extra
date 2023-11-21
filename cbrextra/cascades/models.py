from sqlalchemy import Column, Float, ForeignKey, Integer, String, Sequence
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()

class OrganismModel(Base):
    __tablename__ = "organisms"
    tax_id = Column(Integer, primary_key=True)
    name = Column(String)
    cascade_organisms = relationship('CascadeOrganism', back_populates='organism')

class CascadeStepOrganismModel(Base):
    __tablename__ = "cascade_step_organism"
    id = Column(Integer, Sequence("cascade_organism_id"), primary_key=True)
    identity = Column(Float)
    sequence_match = Column(String)

    organism_id = Column(Integer, ForeignKey("organisms.tax_id"))
    organism = relationship("OrganismModel", back_populates='cascade_organisms')
    
    step_id = Column(Integer, ForeignKey('cascade_steps.id'))
    step = relationship('CascadeStepModel', back_populates='step_id')

class CascadeStepSequenceModel(Base):
    __tablename__ = "cascade_step_sequences"
    id = Column(Integer, Sequence("cascade_sequence_id"), primary_key=True)

    # Sequence
    seq_id = Column(String)
    sequence = Column(String)

    # Cascade Step
    step_id = Column(Integer, ForeignKey('cascade_steps.id'))
    step = relationship('CascadeStepModel', back_populates='step_id')

class CascadeStepModel(Base):
    __tablename__ = "cascade_steps"
    id = Column(Integer, primary_key=True)
    step_name = Column(String)
    step_sequences = relationship('CascadeStepSequenceModel', back_populates='step')
    step_organisms = relationship('CascadeStepOrganismModel', back_populates='step')