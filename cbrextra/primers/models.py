from sqlalchemy import Column, Float, ForeignKey, Integer, String, Sequence, Text, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from typing import Any, Dict, Iterable

from .data import DesignPrimersResults, PrimerResult

Base = declarative_base()

class CodonMappingModel(Base):
    __tablename__ = 'codon_mappings'
    __table_args__ = (UniqueConstraint("plasmid_id", "amino_acid", name="unambigous_coding_constraint"),)
    id = Column(Integer, Sequence('codon_mapping_id'), primary_key=True)
    plasmid_id = Column(Integer, ForeignKey("plasmids.id"))
    plasmid = relationship("PlasmidModel", back_populates="mappings")
    codon = Column(String(3))
    amino_acid = Column(String(3))

    @staticmethod
    def from_results(plasmid : 'PlasmidModel', results : DesignPrimersResults) -> 'Iterable[CodonMappingModel]':

        for aa, codon in results.codon_mappings.items():
            yield CodonMappingModel(
                plasmid=plasmid,
                codon=codon,
                amino_acid=aa
            )

class PlasmidModel(Base):
    __tablename__ = 'plasmids'
    id = Column(Integer, Sequence('plasmid_id'), primary_key=True)
    sequence = Column(Text)
    primers = relationship("PrimerModel", back_populates="plasmid")
    mappings = relationship("CodonMappingModel", back_populates="plasmid")

    @staticmethod
    def from_results(results : DesignPrimersResults) -> 'PlasmidModel':
        return PlasmidModel(
            sequence = results.plasmid
        )

    def get_codon_mappings(self) -> Dict[str, str]:
        return dict(
            (mapping.amino_acid, mapping.codon)
            for mapping in self.mappings
        )

class PrimerModel(Base):
    __tablename__ = 'primers'
    id = Column(Integer, Sequence('primier_id'), primary_key=True)
    primer_left = Column(String(50))
    primer_right = Column(String(50))
    tm_left = Column(Float)
    tm_right = Column(Float)
    inner_seq = Column(String(50))
    tm_all = Column(Float)
    position = Column(Integer)
    amino_acid = Column(String(3))
    plasmid_id = Column(Integer, ForeignKey("plasmids.id"))
    plasmid = relationship("PlasmidModel", back_populates="primers")

    @staticmethod
    def from_result(plasmid : PlasmidModel, result: PrimerResult) -> 'PrimerModel':
        return PrimerModel(
            primer_left = result.left_primer,
            primer_right = result.right_primer,
            tm_left = result.tm_left,
            tm_right = result.tm_right,
            inner_seq = result.inner_seq,
            tm_all = result.tm_all,
            position = result.position,
            amino_acid = result.amino_acid,
            plasmid = plasmid
        )

    def to_result(self : Any) -> PrimerResult:
        return PrimerResult(
            left_primer=self.primer_left,
            right_primer=self.primer_right,
            tm_left=self.tm_left,
            tm_right=self.tm_right,
            tm_all=self.tm_all,
            position=self.position,
            amino_acid=self.amino_acid,
            inner_seq=self.inner_seq
        )

    @staticmethod
    def to_results(primers_seq : 'Iterable[PrimerModel]') -> DesignPrimersResults:
        primers = list(primers_seq)
        plasmid = primers[0].plasmid
        plasmid_seq = plasmid.sequence
        mappings = plasmid.get_codon_mappings()

        def convert(model : 'PrimerModel') -> PrimerResult:
            assert model.plasmid.sequence == plasmid_seq, "Corrupt database! Primers have more than one plasmid"
            return model.to_result()

        return DesignPrimersResults(
            plasmid=plasmid_seq,
            codon_mappings=mappings,
            primers = [convert(primer) for primer in primers]
        )
