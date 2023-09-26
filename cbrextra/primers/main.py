"""
usage:
    cbrtools primers design
"""

from docopt import docopt

from ..core.module import Context, Module, Result

class PrimerModule(Module):

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if options.get('design'):
            return self.success

        return self.not_requested