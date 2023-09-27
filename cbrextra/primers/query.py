import sys
from typing import TextIO

from .store import Store

def query_best_primers(
        database: str,
        tm_all : float,
        tm_all_weight : float,
        tm_primers : float,
        tm_primers_weight : float,
        tm_delta_weight : float,
        out_stream : TextIO = sys.stdout
    ):
        results = Store.query_best_primers(
            database,
            tm_all,
            tm_all_weight,
            tm_primers,
            tm_primers_weight,
            tm_delta_weight
        )

        out_stream.write(results.model_dump_json())