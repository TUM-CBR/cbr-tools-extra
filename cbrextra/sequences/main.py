"""Usage:
    cbrtools sequences scan --db-file=<db-file> [-- <scan-dir>...]
    cbrtools sequences tblastn --db-file=<db-file>

Environmental variables:
    CBR_MAKEBLAST_DB    The location of the 'makeblastdb' command.
    CBR_TBLASTN         The location of the 'tblastn' command.
"""

from cbrextra.core.module import Context, Result
from docopt import docopt
from typing import Optional, Sequence, TextIO, Union
from os import environ
import sys

from ..core.module import Module
from .data import BlastEnv, SequencesContext
from .procedures import run_query_tblastn, run_scan

K_MAKEBLAST_DB = "CBR_MAKEBLAST_DB"
K_TBLASTN = "CBR_TBLASTN"

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

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if options.get("scan"):
            return self.scan(options['--db-file'], options.get('<scan-dir>'))
        elif options.get("tblastn"):
            return self.query(options['--db-file'])
        
        return Result.not_requested()
    
module = SequencesModule()