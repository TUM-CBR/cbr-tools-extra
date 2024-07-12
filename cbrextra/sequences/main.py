"""Usage:
    cbrtools sequences scan --db-file=<db-file> [-- <scan-dir>...]
    cbrtools sequences tblastn --db-file=<db-file>
    cbrtools sequences errors --db-file=<db-file>
    cbrtools sequences interactive [--email=<email>] [--api-key=<api_key>]

Environmental variables:
    CBR_MAKEBLAST_DB    The location of the 'makeblastdb' command.
    CBR_TBLASTN         The location of the 'tblastn' command.
"""

from Bio import Entrez
from docopt import docopt
from typing import Optional, Sequence, TextIO, Union
from os import environ
import sys

from ..core.interactive import InteractiveSpec, OnMessageArgs
from ..core.module import Context, Module, Result
from .data import BlastEnv, InteractiveInput, InteractiveOutput, SequencesContext
from .procedures import run_query_errors, run_query_tblastn, run_scan
from .search import Search, SearchArgs

K_MAKEBLAST_DB = "CBR_MAKEBLAST_DB"
K_TBLASTN = "CBR_TBLASTN"

class SequencesInteractive(InteractiveSpec[InteractiveInput, InteractiveOutput]):

    def __init__(self):
        self.__ncbi_search = Search()

    def __search(self, search_args: SearchArgs):
        return self.__ncbi_search.run(search_args)

    def on_message(self, args: OnMessageArgs[InteractiveInput, InteractiveOutput]) -> None:

        message = args.message

        if message.search is not None:

            for result in self.__search(message.search):
                args.send(
                    InteractiveOutput(
                        search_result=result
                    )
                )
                return

        raise ValueError("Message does not contain a valid command")

class SequencesModule(Module):

    def __blast_context_from_env(self) -> BlastEnv:

        return BlastEnv(
            makeblastdb=environ.get(K_MAKEBLAST_DB, "makeblastdb"),
            tblastn=environ.get(K_TBLASTN, "tblastn")
        )

    def scan(self, db_file: str, scan_dirs: Optional[Union[str, Sequence[str]]]):

        if scan_dirs is None:
            scan_dirs = []
        elif isinstance(scan_dirs, str):
            scan_dirs = [scan_dirs]

        run_scan(
            SequencesContext(
                self.__blast_context_from_env(),
                db_file=db_file
            ),
            scan_dirs
        )

        return Result.success()
    
    def query(
        self,
        db_file: str,
        in_stream: TextIO = sys.stdin,
        out_stream: TextIO = sys.stdout
    ):
        run_query_tblastn(
            SequencesContext(
                self.__blast_context_from_env(),
                db_file=db_file
            ),
            query=in_stream,
            result_stream=out_stream
        )

        return Result.success()

    def setup_entrez(
        self,
        email: Optional[str] = None,
        api_key: Optional[str] = None
    ):
        Entrez.email = email
        Entrez.api_key = api_key

    def interactive_mode(
        self,
        email: Optional[str],
        api_key: Optional[str],
        in_stream: TextIO = sys.stdin,
        out_stream: TextIO = sys.stdout
    ) -> Result:

        self.setup_entrez(email, api_key)
        search.run(
            email,
            api_key,
            in_stream = in_stream,
            out_stream = out_stream
        )

        return Result.success()
    
    def errors(
        self,
        db_file: str,
        out_stream: TextIO = sys.stdout
    ):
        run_query_errors(
            SequencesContext(
                self.__blast_context_from_env(),
                db_file=db_file
            ),
            result_stream=out_stream
        )

        return Result.success()

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if options.get("interactive"):
            return self.interactive_mode(
                email=options.get("--email"),
                api_key=options.get("--api-key")
            )
        if options.get("scan"):
            return self.scan(options['--db-file'], options.get('<scan-dir>'))
        elif options.get("tblastn"):
            return self.query(options['--db-file'])
        elif options.get("errors"):
            return self.errors(options['--db-file'])
        
        return Result.not_requested()
    
module = SequencesModule()
