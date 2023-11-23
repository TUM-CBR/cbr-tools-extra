from operator import and_, or_
from Bio.SeqRecord import SeqRecord
from os import path
from typing import Iterable, Tuple
from sqlalchemy import Select, create_engine, select, join
from sqlalchemy.orm import Session, sessionmaker

from cbrextra.cascades.find_organisms import default_exclude

from .data import *
from .models import Base, CascadeStepModel, CascadeStepOrganismModel, CascadeStepSequenceModel, OrganismModel

class Store:

    class StoreApi:

        def __init__(self, session: Session):
            self.__session = session

        def __enter__(self):
            self.__session.__enter__()
            return self

        def __exit__(self, *args, **kwargs):
            self.__session.commit()
            self.__session.__exit__(*args, **kwargs)

        def query_organisms_for_step(self, step : CascadeStepModel) -> Iterable[OrganismModel]:
            table = join(
                OrganismModel,
                CascadeStepOrganismModel,
                OrganismModel.tax_id == CascadeStepOrganismModel.organism_id
            )
            query = select(OrganismModel) \
                .where(CascadeStepOrganismModel.step_id == step.id) \
                .select_from(table)

            return [
                organism
                for (organism,) in self.__session.execute(query)
            ]

        def query_sequences_for_step(self, step : CascadeStepModel) -> Iterable[CascadeStepSequenceModel]:
            query = select(CascadeStepSequenceModel).where(CascadeStepSequenceModel.step_id == step.id)

            return [
                seq
                for (seq,) in self.__session.execute(query)
            ]

        def load_cascade_args(self) -> Iterable[FindOrganismsArgs]:

            for (step,) in self.__session.execute(select(CascadeStepModel)):

                previous_organisms = [
                    f"{organism.name} (taxid:{organism.tax_id})"
                    for organism in self.query_organisms_for_step(step)
                ]
                excluded_organisms = default_exclude + previous_organisms

                sequences = [
                    CascadeSequence(seq_id = seq.seq_id, seq = seq.sequence)
                    for seq in self.query_sequences_for_step(step)
                ]

                yield FindOrganismsArgs(
                    step = CascadeStep(
                        step_id = step.id,
                        step_name = step.step_name,
                        sequences = sequences
                    ),
                    excluded_organisms=excluded_organisms
                )


            query = select(CascadeStepModel)
            return self.__session.execute(query)

        def save_step_organisms(self, step : CascadeStepModel, step_organisms: Iterable[CascadeStepOrganism]):

            self.__session.add_all(
                CascadeStepOrganismModel(
                    sequence_match = step_organism.sequence_match,
                    step = step,
                    organism = self.get_or_create_organism(step_organism.organism),
                    identity = step_organism.identity
                )
                for step_organism in step_organisms
            )

        def get_or_create_organism(self, organism: Organism) -> OrganismModel:
            existing_organism = self.__session.execute(
                select(OrganismModel).where(OrganismModel.tax_id == organism.tax_id)
            ).first()

            if existing_organism is not None:
                return existing_organism[0]

            new_model = OrganismModel(
                tax_id = organism.tax_id,
                name = organism.name
            )

            self.__session.add(new_model)
            return new_model
            

        def save_results(self, results: Iterable[CascadeStepResult]) -> None:
            for result in results:
                step_result = self.__session.execute(select(CascadeStepModel) \
                    .where(CascadeStepModel.id == result.step.step_id)) \
                    .first()

                assert step_result is not None
                (step,) = step_result

                self.save_step_organisms(step, result.organisms)

        def create_step(self, value: dict) -> CascadeStepModel:
            model = CascadeStepModel.from_dict(value)
            self.__session.add(model)
            return model

        def create_step_seq(self, step : CascadeStepModel, seq: SeqRecord):
            model = CascadeStepSequenceModel(
                seq_id = seq.id,
                sequence = str(seq.seq),
                step = step
            )
            self.__session.add(model)
            return model

        def get_steps(self, step_ids: List[int]) -> List[CascadeStepModel]:
            query = select(CascadeStepModel).where(
                CascadeStepModel.id.in_(step_ids)
            )

            results = list(self.__session.execute(query))
            for (result,) in results:
                if result.id not in step_ids:
                    raise ValueError(f"The step {step_ids} is not in this database.")

            return [result for (result,) in results]

        def get_organisms_missing_step_gene(
            self,
            step : CascadeStepModel,
            max_identity_treshold : float
        ) -> List[Tuple[OrganismModel, float]]:

            if not (0 <= max_identity_treshold <= 1):
                raise ValueError(f"The 'max_identity_treshold' should be between 0 and 1. Got {max_identity_treshold}")

            table = join(
                OrganismModel,
                CascadeStepOrganismModel,
                OrganismModel.tax_id == CascadeStepOrganismModel.organism_id
            )

            # we find all the organisms that have genes similar
            # to the sequence of the enzymes of the provided step
            step_organisms_ids = \
                select(CascadeStepOrganismModel.organism_id) \
                .where(CascadeStepOrganismModel.step_id == step.id)

            query_low_identity = \
                select(OrganismModel, CascadeStepOrganismModel).select_from(table) \
                .where(
                    # The identity between genes of the organism is low
                    # aganist the gene sequence of the enzyme for this step
                    and_ (
                        CascadeStepOrganismModel.step_id == step.id,
                        CascadeStepOrganismModel.identity <= max_identity_treshold
                    )
                )

            query_missing = \
                select(OrganismModel) \
                .where(
                    # Organisms that don't have any genes similar
                    # to the gene of the enzymes in the step
                    CascadeStepOrganismModel.organism_id.not_in(step_organisms_ids)
                )

            low_identity = [
                (organism, identity)
                for (organism, identity) in self.__session.execute(query_low_identity)
            ]

            missing = [
                (organism, 0)
                for (organism,) in self.__session.execute(query_missing)
            ]

            return low_identity + missing

        def query_step_identities(
            self,
            step_id : int,
            organisms: Iterable[OrganismModel]
        ) -> List[CascadeStepOrganismModel]:

            ids = (organism.tax_id for organism in organisms)
            query = Select(CascadeStepOrganismModel) \
                .where(
                    and_(
                        CascadeStepOrganismModel.id.in_(ids),
                        CascadeStepSequenceModel.step_id == step_id
                    )
                )

            return [
                result
                for (result,) in self.__session.execute(query)
            ]

    def __init__(self, db_file: str):
        self.__db_file = db_file
        self.__engine = create_engine(f"sqlite:///{self.__db_file}")
        self.__init_db()
        self.__session_factory = sessionmaker(bind = self.__engine)

    def __init_db(self):
        if path.exists(self.__db_file):
            return

        Base.metadata.create_all(self.__engine)

    def session(self):
        return Store.StoreApi(self.__session_factory())