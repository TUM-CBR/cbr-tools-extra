"""Usage:
    cbrtools primers design <output>
    cbrtools primers query [--tm_total=<tm_total_target>,<tm_total_weight>] [--tm_delta=<tm_delta_target>,<tm_delta_weight>] [--tm_primers=<tm_primers_target>,<tm_priemrs_weight>] <primers_db>

Options:
    --tm_total=<tm_total_target>,<tm_total_weight>                   # Tuple specifying the taget tm of the whole amplicon and the
                                                                     # weight given to that when scoring.
    --tm_delta_weight=<tm_delta_weight>                              # Weight given to the tm difference between the left and
                                                                     # right primer when scoring.
    --tm_primers=<tm_primers_target>,<tm_priemrs_weight>             # Tuple specifying the target tm each of the primers should
                                                                     # have and the weight given to that parameter when scoring.
    <primers_db>                                                     # Primers database to query the values

"""
from docopt import docopt
import json
import os
import sys
from typing import Any, Dict, Tuple

from ..core.module import Context, Module, Result
from .design import DesignPrimersArgs, Operations
from .store import Store
from . import query

class PrimerModule(Module):

    def __design_primers(self, options : Dict[str, Any]) -> Result:

        out_file = options.get("<output>")

        if out_file is None:
            return self.invalid_arguments({"--output": "Missing results output file."})

        if os.path.exists(out_file):
            os.remove(out_file)

        input_args = json.load(sys.stdin)
        args = DesignPrimersArgs(**input_args)
        primers = Operations.design_primers(args)
        Store.save_design_results(primers, out_file)

        return self.success

    def __read_value_and_weight(self, option: str, options: Dict[str, Any]) -> Tuple[float, float]:
        value = options.get(option)
        if value is None:
            return (0,0)

        values = value.split(',')
        error = ValueError(f"The option '{option}' must have the format '#.#,#.#'")

        if len(values) != 2:
            raise error

        try:
            return (float(values[0]), float(values[1]))
        except ValueError:
            raise error

    def __query_primers(self, options : Dict[str, Any]) -> Result:

        tm_total, tm_total_weight = self.__read_value_and_weight("--tm_total", options)
        tm_primers, tm_primers_weight = self.__read_value_and_weight("--tm_primers", options)

        try:
            tm_delta_weight = float(options.get("--tm_delta_weight") or "1")
        except ValueError:
            raise ValueError("The option '--tm_delta_weight' must be a real number.")

        dataset = options['<primers_db>']

        query.query_best_primers(
            dataset,
            tm_total,
            tm_total_weight,
            tm_primers,
            tm_primers_weight,
            tm_delta_weight
        )

        return self.success

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if options.get('design'):
            return self.__design_primers(options)
        if options.get('query'):
            return self.__query_primers(options)

        return self.not_requested

module = PrimerModule()