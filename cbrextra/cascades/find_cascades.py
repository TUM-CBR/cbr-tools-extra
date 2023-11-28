import asyncio
from Bio import SeqIO
from concurrent.futures import Executor, ThreadPoolExecutor
from typing import cast, List

from .find_organisms import find_organisms
from .store import Store

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

MIN_RESULTS_UNTIL_FAILURE = 10

async def add_missing(
    step_id: int,
    store : Store,
    executor : Executor,
    throw_exceptions = False
    ):

    with store.session() as session:
        missing = list(session.load_missing_args([step_id]))

    assert len(missing) < 2
    if len(missing) == 0:
        return

    missing_arg = missing[0]
    missing_result = None

    try:
        missing_result = await find_organisms(missing_arg, executor)
    except Exception as e:
        if throw_exceptions:
            raise e

    if missing_result is not None:
        with store.session() as session:
            session.save_results([missing_result])

async def build_cascades_for_step(
    step_id : int,
    store : Store,
    target_identity : float,
    executor : Executor
):

    with store.session() as session:
        [arg] = session.load_cascade_args([step_id])

    current_identity = 1
    num_results = arg.num_results
    retry_delay = 1

    while(current_identity > target_identity):

        try:
            results = await find_organisms(
                arg._replace(num_results = num_results),
                executor
            )
        except TypeError as e:
            raise e
        except Exception as e:
            num_results = num_results / 2
            retry_delay *= 2
            if num_results < MIN_RESULTS_UNTIL_FAILURE:
                raise e
            await asyncio.sleep(retry_delay)
            continue

        current_identity = \
            float(0) if len(results.organisms) == 0 else \
            min(*[organism.identity for organism in results.organisms])

        with store.session() as session:
            session.save_results([results])
            [arg] = session.load_cascade_args([step_id])

        await add_missing(step_id, store, executor)

def build_cascades_db(
    results_db : str,
    target_identity : float
):

    store = Store(results_db)
    
    with store.session() as session:
        steps = [cast(int, step.id) for step in session.get_steps()]

    with ThreadPoolExecutor(max_workers=4) as executor:
        operations = [
            build_cascades_for_step(cast(int, step), store, target_identity, executor)
            for step in steps
        ]

        # If a single operation fails, don't cancel the
        # rest. Collect the results and exceptions and
        # raise afterwards if necessary
        results = asyncio.gather(
            *operations,
            return_exceptions=True,
        )

        for result in results.result():
            if isinstance(result, BaseException):
                raise result

        # Query missing one last time with all the results
        # that have been gathered
        asyncio.gather(*[add_missing(step, store, executor) for step in steps]).result()