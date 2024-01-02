"""Usage:
    cbrtools kinetics eval
"""

from docopt import docopt
import json
import sys
from typing import Any, Dict

from ..core.module import Context, Module, Result
from .data import EvalArgs

class EnzymeKinetics(Module):

    def __eval_model(self, options: Dict[str, Any]) -> Result:
        input_args = json.load(sys.stdin)
        args = EvalArgs(**input_args)

        

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if options.get('eval'):
            return self.__eval_model(options)