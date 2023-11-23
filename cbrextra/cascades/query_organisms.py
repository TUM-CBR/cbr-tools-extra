from .data import *
from .store import OrganismModel, Store, Tuple
from typing import cast

def find_cascades(
    cascade_args: QueryCascadeArgs,
    db_file: str
) -> QueryCascadeResult:
    """Search in a database of cascades and organisms for organisms of interest.
    In the arguments, a maximum identity treshold is provided which sets the maximum
    limit at which two genes are to be considered the same. A list of steps is also
    provided and for each step one must specify wether:
        * keep: indicates that one is only interested in organisms that have genes
                that encode any of the enzymes needed for this step.
        * replace: indicates one is only interested in organisms that don't have genes
                   that code for the any of the enzymes needed for the step
        * any: 'keep' or 'replace'
    
    This function will query all the organisms mathcing all of those steps and filter
    them out according to the criteria specified above.
    """

    store = Store(db_file)
    steps_specs = cascade_args.steps
    by_step = dict((arg.step_id, arg) for arg in steps_specs)

    with store.session() as session:
        steps = session.get_steps(list(by_step.keys()))
        step_identities = session.get_organisms_and_steps(steps)

    organisms = dict(
        (organism.tax_id, organism)
        for (_, organism) in step_identities
    )

    all_steps : Dict[int, Dict[int, float]] = {}

    # todo: This is very inefficient as we are not using the db manager at all
    # to do the filtering. Currently, this is not a problem as the results come
    # from NCBI, so we are anyways limited by their throttling regarding how
    # big these databases can become.
    # the commit c1480b40283f5b7caaf19b116bcab2c1a336e8a6 contains a start
    # for that implementation
    for step_organism, organism in step_identities:
        step_id = cast(int, step_organism.step_id)
        all_steps[step_id] = step_dict = all_steps[step_id] if step_id in all_steps else {}
        tax_id = cast(int, organism.tax_id)

        step_dict[tax_id] = cast(float, step_organism.identity)

    def construct_steps(organism: OrganismModel):

        result = []
        for step in steps_specs:
            identity = next(
                (
                    identity
                    for identities in [all_steps.get(step.step_id)] if identities is not None
                    for identity in [identities.get(cast(int, organism.tax_id))]
                ),
                None
            )

            # ideally, these two filtering steps should be done by
            # the dbms
            if identity is None and step.policy == QueryStepPolicy.keep:
                return None

            if identity is not None \
                and step.policy == QueryStepPolicy.replace \
                and identity > cascade_args.max_identity_treshold:
                return None
            
            result.append(
                QueryCascadeResultStepEntry(
                    step_id=step.step_id,
                    identity=identity or 0
                )
            )

    result_organisms = [
        QueryCascadeResultOrganismEntry(
            organism = Organism(tax_id=cast(int, organism.tax_id), name = cast(str, organism.name)),
            steps = steps
        )
        for organism in organisms.values()
        for steps in [construct_steps(organism)] if steps is not None
    ]

    return QueryCascadeResult(organisms=result_organisms)