"""CBR Tools Extras

Usage:
    cbrtools.py [--help] [--version] <command> [<args>...]
  
Options:
  -h --help     Show this screen.
  --version     Show version.

The most commonly used commands are:
   energy        Perform or read energy calculations
"""
from sys import argv, exit, stderr
from typing import Optional
from docopt import docopt

from cbrextra.core.module import Context, Result
from cbrextra.cascades.main import module as cascade
from cbrextra.energy.main import module as energy
from cbrextra.primers.main import module as primers

BAD_OPERATION = [
    "The request:",
    "",
    " ".join(argv),
    "",
    "was not understood by 'cbrtools'. The supported operations are:",
    "",
    __doc__
]

def main():
    arguments = docopt(__doc__, version='CBR Tools Extra 1.0', options_first=True)
    context = Context(cmdargs=arguments)
    result : Optional[Result] = None
    command = arguments.get("<command>")
    
    if command == 'energy':
        result = energy.main(context)
    elif command == 'primers':
        result = primers.main(context)
    elif command == 'cascades':
        result = cascade.main(context)
    #elif command == "kinetics":
    #    from cbrextra.kinetics.main import module as kinetics
    #    result = kinetics.main(context)
    elif command == "cavities":
        from cbrextra.cavities.main import module as cavities
        result = cavities.main(context)
    elif command == "coevolution":
        from cbrextra.coevolution.main import module as coevolution
        result = coevolution.main(context)
    elif command == "sequences":
        from cbrextra.sequences.main import module as sequences
        result = sequences.main(context)

    if result is None or not result.is_success:
        print(docopt(__doc__), file=stderr)
        exit(1)