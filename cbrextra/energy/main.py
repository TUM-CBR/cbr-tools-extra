"""
usage:
    cbrtools energy parseresults <directory> [--out-dir=<directory>]

Generic options
    -h, --help

Specific energy options:
    --out-dir=<directory>   Directory to store the parsed results
"""

from docopt import docopt

from ..core.module import Context, Module, Result

from .postprocessing import ParseEnergyRunsOperation

class EnergyModule(Module):

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        parseresults = options.get('parseresults')
        directory = options.get("<directory>")
        o_directory = options.get("--out-dir")
        if parseresults and directory:
            ParseEnergyRunsOperation(context).to_csv(
                directory,
                o_directory
            )

            return self.success

        else:
            return self.not_requested

module = EnergyModule()