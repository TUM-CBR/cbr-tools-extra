"""Usage:
    cbrtools cascades query <database> <max-identity-treshold> [<step-id>,<policy>...]
    cbrtools cascades create <target-identity> <email> <database> <fasta-file> <spec-file>
"""
from Bio import Entrez
from docopt import docopt
import json
from os import path
import re
import sys
from typing import List, TextIO

from ..core.module import Context, Module, Result

from .data import QueryCascadeArgs, QueryCascadeStep, QueryStepPolicy
from .find_cascades import initialize_cascade_database, build_cascades_db
from .query_organisms import find_cascades


email_pattern = re.compile(r"^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$")

def is_email(email):
    return email_pattern.match(email) is not None

INVALID_EMAIL_MSG = """You must provide a valid email to use this feature. This email
is provided to NCBI on any request. As this program queries the NCBI BLAST service
quite aggressively, NCBI will notify you by email if you are using the resources
too excesively.
"""

class CascadesModule(Module):

    def __create_cascades(
        self,
        target_identity_str : str,
        email : str,
        db_file : str,
        fasta_file : str,
        spec_file_name : str
    ):

        if not is_email(email):
            raise ValueError(INVALID_EMAIL_MSG)

        Entrez.email = email
        target_identity = float(target_identity_str)

        if not (0 < target_identity <= 1):
            raise ValueError(f"Identity must be a number between 0 and 1, got {target_identity}")

        if path.exists(db_file):
            raise ValueError(
            f"Cowardly refusing to overwrite database at '{db_file}'"
        )
    
        with open(spec_file_name, 'w') as spec_file:
            spec = json.load(spec_file)


        initialize_cascade_database(
            db_file,
            fasta_file,
            spec
        )

        build_cascades_db(
            db_file,
            target_identity
        )

        return Result.success()

    def __query_cascades(
        self,
        db_file : str,
        max_identity_treshold : str, 
        steps : List[str],
        out_stream : TextIO = sys.stdout
    ):
        query_cascade_args = [
            QueryCascadeStep(
                step_id=int(step_id),
                policy=QueryStepPolicy.read(policy)
            )
            for step_id,policy in (step.split() for step in steps)
        ]

        result = find_cascades(
            QueryCascadeArgs(
                steps = query_cascade_args,
                max_identity_treshold = float(max_identity_treshold)
            ),
            db_file=db_file
        )

        out_stream.write(result.model_dump_json())

        return Result.success()

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if not options.get('cascades'):
            return Result.not_requested()

        if options.get('query'):
            db_file = options['<database>']
            max_identity_treshold = options['<max-identity-treshold>']
            steps = options['<step-id>,<policy>']

            return self.__query_cascades(
                db_file,
                max_identity_treshold,
                steps
            )
        elif options.get('create'):
            db_file = options['<database>']
            target_identity = options['<target-identity>']
            fasta_file = options['<fasta-file>']
            spec_file = options['<spec-file>']
            email = options['<email>']

            return self.__create_cascades(
                target_identity,
                email,
                db_file,
                fasta_file,
                spec_file
            )

module = CascadesModule()