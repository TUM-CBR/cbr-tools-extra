from itertools import islice
from sqlalchemy import create_engine, func, select
from sqlalchemy.orm import sessionmaker
from typing import List, NamedTuple

from .data import DesignPrimersResults
from .models import Base, CodonMappingModel, PlasmidModel, PrimerModel

class PrimersQueryResult(NamedTuple):
    primers : List[PrimerModel]
    plasmid : PlasmidModel

class Store:

    def __init__(self, db_file: str):
        self.__db_file = db_file

    def __create_engine(self):
        return create_engine(f"sqlite:///{self.__db_file}")

    def __query_best_primers(
        self,
        tm_all : float,
        tm_all_weight : float,
        tm_primers : float,
        tm_primers_weight : float,
        tm_delta_weight : float
    ):

        tm_primers_weight = tm_primers_weight / 2

        Session = sessionmaker(bind=self.__create_engine())

        query = select(
            PrimerModel,
            func.min(
                tm_all_weight*func.abs(PrimerModel.tm_all - tm_all) \
                + func.abs(PrimerModel.primer_left - PrimerModel.primer_right)*tm_delta_weight \
                + func.abs(PrimerModel.primer_left - tm_primers)*tm_primers_weight \
                + func.abs(PrimerModel.primer_right - tm_primers)*tm_primers_weight
            )
        ).group_by('amino_acid', 'position')

        with Session() as session:
            return PrimerModel.to_results(
                primer
                for primer,_ in session.execute(query)
            )

    def __save_design_results(
        self,
        primer_results: DesignPrimersResults,
    ):
        db_engine = self.__create_engine()
        Base.metadata.create_all(db_engine)

        Session = sessionmaker(bind=db_engine)

        with Session() as session:
            plasmid = PlasmidModel.from_results(primer_results)
            codons = CodonMappingModel.from_results(plasmid, primer_results)
            session.add(plasmid)
            session.add_all(codons)
            session.commit()

        all_primers = (
            PrimerModel.from_result(plasmid, primer)
            for primer in primer_results.primers
        )

        while len(chunk := list(islice(all_primers, 10000))) > 0:
            with Session() as session:
                session.add_all(chunk)
                session.commit()

    @staticmethod
    def save_design_results(
        primer_results: DesignPrimersResults,
        out_file: str
    ):
        Store(out_file).__save_design_results(primer_results)

    @staticmethod
    def query_best_primers(
        db: str,
        tm_all : float,
        tm_all_weight : float,
        tm_primers : float,
        tm_primers_weight : float,
        tm_delta_weight : float
    ):
        return Store(db).__query_best_primers(
            tm_all,
            tm_all_weight,
            tm_primers,
            tm_primers_weight,
            tm_delta_weight
        )
