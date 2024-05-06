import os
import subprocess

from .data import BlastEnv

class Blast:

    def __init__(
        self,
        env: BlastEnv
    ):
        
        self.__env = env

    @property
    def __flags(self) -> int:

        if os.name == "nt":
            flags = subprocess.CREATE_NO_WINDOW
        else:
            flags = 0

        return flags

    def update_dbs(
        self,
        database: str,
        *fasta_files: str
    ):
        
        for fasta in fasta_files:

            with subprocess.Popen(
                [
                    self.__env.makeblastdb,
                    "-in",
                    fasta,
                    "-dbtype", "nucl",
                    "-out",
                    database
                ],
                text=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                creationflags=self.__flags
            ) as blast:
                blast.wait()