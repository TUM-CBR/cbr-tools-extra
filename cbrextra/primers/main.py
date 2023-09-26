"""
usage:
cbrtools primers design <output>

"""
from docopt import docopt
import json
import os
import sys
from typing import Any, Dict

from ..core.module import Context, Module, Result
from .design import DesignPrimersArgs, Operations
from .store import Store

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

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if options.get('design'):
            return self.__design_primers(options)

        return self.not_requested

module = PrimerModule()