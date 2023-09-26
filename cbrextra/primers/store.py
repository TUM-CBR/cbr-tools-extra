from itertools import islice
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from .data import DesignPrimersResults
from .models import Base, PrimerModel

class Store:

    def __save_design_results(
        self,
        primer_results: DesignPrimersResults,
        out_file: str
    ):

        db_engine = create_engine(f"sqlite:///{out_file}")
        Base.metadata.create_all(db_engine)

        Session = sessionmaker(bind=db_engine)
        all_primers = (
            PrimerModel.from_result(group_key, codon, result)
            for group_key, group in enumerate(primer_results)
            for codon, results in group.items()
            for result in results
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
        Store().__save_design_results(primer_results, out_file)
