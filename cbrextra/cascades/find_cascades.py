import asyncio
import time
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

) -> List[CascadeStepResult | BaseException]:

    with ThreadPoolExecutor(max_workers=4) as executor:
        coro : Awaitable[List[CascadeStepResult | BaseException]] = asyncio.gather(*
                [
                    find_organisms(arg, executor)
                    for arg in args
                ],
                return_exceptions=True
            )
        results = asyncio.get_event_loop().run_until_complete(coro)

    return results

MIN_RESULTS_UNTIL_FAILURE = 10

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

    # In case of a failure, we wait some time to
    # avoid over-using the NCIB resources
    retry_delay = 1

    while(len(args) > 0):

        results_with_exn = find_cascades_step(args)

        failed = {}
        success = []

        for i, outcome in enumerate(results_with_exn):

            # If we encounter a failre, we retry the query with
            # asking for less results in return
            if isinstance(outcome, BaseException):
                arg = args[i]
                num_results = int(arg.num_results / 2)

                if num_results < MIN_RESULTS_UNTIL_FAILURE:
                    raise outcome
                failed[arg.step.step_id] = arg._replace(num_results = num_results)

            else:

                identity = \
                    float(0) if len(outcome.organisms) == 0 else \
                    min(*[organism.identity for organism in outcome.organisms])

                identities[outcome.step.step_id] = identity
                success.append(outcome)

        with store.session() as session:
            session.save_results(success)
            args = [
                arg
                for arg in session.load_cascade_args()
                if arg.step.step_id not in failed
                    and identities[arg.step.step_id] > target_identity
            ] + list(failed.values())

        # Exponentially increase the delay before retrying
        # if we encountered any failures
        retry_delay *= 1 if len(failed) == 0 else 2
        time.sleep(retry_delay - 1)