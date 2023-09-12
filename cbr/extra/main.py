"""CBR Tools Extras

Usage:
    cbrtools.py [--help] [--version] <command> [<args>...]
  
Options:
  -h --help     Show this screen.
  --version     Show version.

The most commonly used git commands are:
   energy        Perform or read energy calculations
"""

from sys import argv
from typing import Optional
from docopt import docopt

from .core.module import Context, Result
from .energy.main import module as energy

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

    if result is None or not result.is_success:
      print()