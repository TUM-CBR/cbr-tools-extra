from io import StringIO
import os
import subprocess
from typing import TextIO

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

            with StringIO() as error_stream \
                , subprocess.Popen(
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
                result = blast.returncode

                assert blast.stderr
                for text in blast.stderr:
                    error_stream.write(text)

                if result != 0:
                    error_stream.seek(0)
                    raise Exception(f"Blast failed: {error_stream.read()}")

    def query_tblastn(self, database: str, query: TextIO, result: TextIO):

        with StringIO() as error_stream \
            , subprocess.Popen(
            [
                self.__env.tblastn,
                "-db", database,
                "-outfmt", "15"
            ],
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        ) as blast:

            assert blast.stdin
            for text in query:
                blast.stdin.write(text)
            blast.stdin.close()

            assert blast.stdout
            for text in blast.stdout:
                result.write(text)

            assert blast.stderr
            for text in blast.stderr:
                error_stream.write(text)

            blast.wait()
            if blast.returncode != 0:
                error_stream.seek(0)
                raise Exception(f"Blast failed: {error_stream.read()}")