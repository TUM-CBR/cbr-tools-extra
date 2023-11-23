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
    
    This function will then first query the database to identify the organisms that satisfy
    the steps with a 'replace' criteria. 
    """

    store = Store(db_file)
    steps_specs = cascade_args.steps
    by_step = dict((arg.step_id, arg) for arg in steps_specs)

    with store.session() as session:
        steps = session.get_steps(list(by_step.keys()))
        organisms_missing_step = dict(
            (
                step.step_id,
                session.get_organisms_missing_step_gene(
                    step,
                    cascade_args.max_identity_treshold
                )
            )
            for step in steps
        )

        organisms : Dict[int, OrganismModel] = dict(
            (organism.tax_id, organism)
            for (_, results) in organisms_missing_step
            for (organism,_) in results
        )

        missing_steps : List[Tuple[int, Dict[int, float]]] = [
            (
                step_id,
                dict(
                    (organism.tax_id, identity)
                    for organism, identity in results
                )
            )
            for step_id, results in organisms_missing_step
        ]

        other_steps : List[Tuple[int, Dict[int, float]]] = [
            (
                step.step_id,
                dict(
                    (tax_id, identity)
                    for cascade_organism in session.query_step_identities(step.step_id, organisms.values())
                    for tax_id, identity in [(cast(int, cascade_organism.organism_id), cast(float, cascade_organism.identity))]
                )
            )
            for step in steps_specs
            if step.policy != QueryStepPolicy.replace
        ]

    all_steps = dict(missing_steps + other_steps)

    def construct_steps(organism: OrganismModel):

        result = []
        for step in steps_specs:
            identity = all_steps[step.step_id].get(cast(int, organism.tax_id))

            if identity is None and step.policy == QueryStepPolicy.keep:
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