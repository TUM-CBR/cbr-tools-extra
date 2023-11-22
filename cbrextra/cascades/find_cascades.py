import asyncio
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from typing import Awaitable, List
from .store import Iterable, Store

from .data import CascadeStepResult, FindOrganismsArgs
from .find_organisms import find_organisms

def initialize_cascade_database(
    results_db : str,
    fasta_file : str,
    json_steps : List[dict]
):
    store = Store(results_db)
    sequences = dict(
        (seq.id, seq)
        for seq in SeqIO.parse(fasta_file, 'fasta')
    )

    with store.session() as session:

        for json_step in json_steps:
            step = session.create_step(json_step)
            for seq_id in json_step["sequences"]:
                seq = sequences[seq_id]
                session.create_step_seq(step, seq)

def find_cascades_step(
    args : Iterable[FindOrganismsArgs]

) -> List[CascadeStepResult]:

    with ThreadPoolExecutor(max_workers=4) as executor:
        coro : Awaitable[List[CascadeStepResult]] = asyncio.gather(*
                [
                    find_organisms(arg, executor)
                    for arg in args
                ]
            )
        results = asyncio.get_event_loop().run_until_complete(coro)

    return results

def find_cascades(
    results_db : str,
    target_identity : float
):

    store = Store(results_db)
    with store.session() as session:
        args = list(session.load_cascade_args())

    identities = dict(
        (arg.step.step_id, float(1))
        for arg in args
    )

    while(any(identity > target_identity for identity in identities.values())):

        results = find_cascades_step(
            arg
            for arg in args if identities[arg.step.step_id] > target_identity
        )

        for result in results:

            identity = \
                float(0) if len(result.organisms) == 0 else \
                min(*[organism.identity for organism in result.organisms])

            identities[result.step.step_id] = identity

        with store.session() as session:
            session.save_results(results)
            args = list(session.load_cascade_args())