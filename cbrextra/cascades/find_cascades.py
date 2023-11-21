import asyncio
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from typing import Awaitable, List
from .store import Store

from .data import CascadeStepResult, Dict
from .find_organisms import find_organisms

def initialize_cascade_database(
    results_db : str,
    fasta_file : str,
    step_mappings : List[dict]
):
    store = Store(results_db)
    sequences = SeqIO.read(fasta_file, 'fasta')
    seq_to_step = dict(
        (seq_id, step_id)
        for step_id, seqs in step_mappings.items()
        for seq_id in seqs
    )

    with store.session() as session:
        steps = dict(
            (step_id, session.create_step(step_id))
            for step_id in step_mappings.keys()
        )

        


def find_cascades_step(
    results_db : str,
) -> List[CascadeStepResult]:
    store = Store(results_db)

    with store.session() as session:
        args = session.load_cascade_args()

    with ThreadPoolExecutor(max_workers=4) as executor:
        coro : Awaitable[List[CascadeStepResult]] = asyncio.gather(*
                [
                    find_organisms(arg, executor)
                    for arg in args
                ]
            )
        results = asyncio.get_event_loop().run_until_complete(coro)

    with store.session() as session:
        session.save_results(results)

    return results