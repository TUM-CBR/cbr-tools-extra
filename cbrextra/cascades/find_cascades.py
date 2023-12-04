import asyncio
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import sys
from typing import cast, List, TextIO

from .find_organisms import find_organisms
from .store import CascadeStepOrganism, Optional, Store
from .support import RunCascadesContext

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
    context : RunCascadesContext,
    step_id: int,
    throw_exceptions = False
    ):

    store = context.store
    executor = context.executor

    with store.session() as session:
        missing = list(session.load_missing_args([step_id]))

    assert len(missing) < 2
    if len(missing) == 0:
        return

    missing_arg = missing[0]
    missing_result = None
    
    context.write_status({
        'step': step_id,
        'missing_organisms': len(missing_arg.included_organisms or [])
    })

    try:
        missing_result = await find_organisms(context, missing_arg, executor)
    except Exception as e:
        context.write_error(step_id, e)
        if throw_exceptions:
            raise e

    if missing_result is not None:
        with store.session() as session:
            session.save_results([missing_result])

MAX_BLAST_CONCURRENT_REQUESTS = 1

UNRECOVERABLE = (TypeError,)

async def add_missing_organisms_main(
    db_file : str,
    out_stream : TextIO,
    included_steps : Optional[List[int]] = None
):
    store = Store(db_file)

    with ThreadPoolExecutor(max_workers = MAX_BLAST_CONCURRENT_REQUESTS) as executor:
        context = RunCascadesContext(
            store = store,
            executor=executor,
            out_stream=out_stream,
            domain=None
        )
        await add_missing_organisms(context, included_steps)

async def add_missing_organisms(
    context : RunCascadesContext,
    included_steps : Optional[List[int]] = None
):

    store = context.store

    with store.session() as session:
        steps = [
            step_id for step in session.get_steps()
            for step_id in [cast(int, step.id)]
            if included_steps is None or step_id in included_steps
        ]

    async def add_missing_retry(step : int):
        retries = 5
        delay = 8

        while True:
            try:
                return await add_missing(context, step, throw_exceptions = True)
            except Exception as e:

                if isinstance(e, UNRECOVERABLE):
                    raise e

                if retries > 0:
                    retries -= 1
                    delay *= 2
                    await asyncio.sleep(delay)
                else:
                    raise e

    # We don't want exceptions to interrupt other steps,
    # nevertheless, we also want to throw the exception (if any)
    exns = await asyncio.gather(
        *[add_missing_retry(step) for step in steps],
        return_exceptions=True
    )

    for e in exns:
        if isinstance(e, BaseException):
            raise e

async def build_cascades_for_step(
    context : RunCascadesContext,
    step_id : int,
    target_identity : float
):
    store = context.store
    executor = context.executor

    with store.session() as session:
        [arg] = context.with_domains(session.load_cascade_args([step_id]))

    current_identity = 1
    num_results = arg.num_results
    retry_delay = 1

    def ping(organisms : Optional[List[CascadeStepOrganism]] = None):

        new_organisms = [('new_organisms', len(organisms))] if organisms is not None else []
        message = dict(
            [
                ('step', step_id),
                ('current_identity', current_identity)
            ]
            + new_organisms
        )

        context.write_status(message)

    ping()

    while(current_identity > target_identity):

        try:
            results = await find_organisms(
                context,
                arg._replace(num_results = num_results),
                executor
            )
        except Exception as e:

            if isinstance(e, UNRECOVERABLE):
                raise e

            context.write_error(
                step_id,
                e,
                {
                    'num_results': num_results
                }
            )
            num_results = num_results / 2
            retry_delay *= 2
            if num_results < MIN_RESULTS_UNTIL_FAILURE:
                raise e
            await asyncio.sleep(retry_delay)
            continue

        current_identity = \
            float(0) if len(results.organisms) == 0 else \
            min(*[organism.identity for organism in results.organisms])

        ping(results.organisms)

        with store.session() as session:
            session.save_results([results])
            [arg] = context.with_domains(session.load_cascade_args([step_id]))

        await add_missing(context, step_id)

async def build_cascades_db(
    results_db : str,
    target_identity : float,
    out_stream : TextIO = sys.stdout,
    domain : Optional[str] = None
):
    store = Store(results_db)
    
    with store.session() as session:
        steps = [cast(int, step.id) for step in session.get_steps()]

    with ThreadPoolExecutor(max_workers=MAX_BLAST_CONCURRENT_REQUESTS) as executor:

        context = RunCascadesContext(
            store = store,
            executor=executor,
            out_stream=out_stream,
            domain=domain
        )

        operations = [
                build_cascades_for_step(
                    context,
                    cast(int, step),
                    target_identity
                )
                for step in steps
        ]

        # If a single operation fails, don't cancel the
        # rest. Collect the results and exceptions and
        # raise afterwards if necessary
        results = await asyncio.gather(*operations, return_exceptions=True)

        for result in results:
            exn = result
            if exn is not None:
                raise exn

        # Query missing one last time with all the results
        # that have been gathered
        await add_missing_organisms(context)