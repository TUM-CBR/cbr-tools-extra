import asyncio
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import sys
from typing import cast, List, TextIO

from .find_organisms import find_organisms
from .store import Optional, Store
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
        missing = context.with_domains(session.load_missing_args([step_id]))

    context.write_status({
        'step': step_id,
        'missing_organisms': len(missing)
    })

    assert len(missing) < 2
    if len(missing) == 0:
        return

    missing_arg = missing[0]
    missing_result = None

    try:
        missing_result = await find_organisms(context, missing_arg, executor)
    except Exception as e:
        context.write_error(step_id, e)
        if throw_exceptions:
            raise e

    if missing_result is not None:
        with store.session() as session:
            session.save_results([missing_result])

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

    def ping():
        context.write_status({
            'step': step_id,
            'current_identity': current_identity
        })

    ping()

    while(current_identity > target_identity):

        try:
            results = await find_organisms(
                context,
                arg._replace(num_results = num_results),
                executor
            )
        except TypeError as e:
            raise e
        except Exception as e:
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

        ping()

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

    with ThreadPoolExecutor(max_workers=1) as executor:

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
        await asyncio.gather(*[add_missing(context, step) for step in steps])