"""Usage:
    cbrtools kinetics eval
"""

from docopt import docopt
import json
import sys
from typing import Any, Dict, TextIO

from ..core.module import Context, Module, Result
from .compute import eval_model
from .data import EvalArgs

class EnzymeKinetics(Module):

    def eval_model(
        self,
        options: Dict[str, Any],
        input_stream: TextIO = sys.stdin,
        output_stream: TextIO = sys.stdout
    ) -> Result:
        input_args = json.load(input_stream)
        args = EvalArgs(**input_args)
        result = eval_model(args)
        output_stream.writelines([
            result.model_dump_json()
        ])
        return Result.success()

    def interactive_mode(
        self,
        options : Dict[str, Any],
        input_stream: TextIO = sys.stdin,
        output_stream: TextIO = sys.stdout
    ) -> Result:


    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if options.get('eval'):
            return self.eval_model(options)

        if options.get('interactive'):
            return self.interactive_mode(options)

        return Result.not_requested()

module = EnzymeKinetics()