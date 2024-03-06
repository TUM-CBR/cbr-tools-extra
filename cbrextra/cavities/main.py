"""Usage:
    cbrtools cavities interactive --input-points=<input_points_json>
"""

from docopt import docopt
from typing import Any, Dict

from cbrextra.core.interactive.shared import OnMessageArgs

from ..core.interactive import *
from ..core.module import Context, Module, Result
from .data import *

class CavitiesInstance:
    pass

class CavititesInteractive(InteractiveSpec[InteractiveInput, InteractiveOutput]):

    def __init__(
        self,
        points: List[Points]
    ):
        pass

    def on_message(self, args: OnMessageArgs[InteractiveInput, InteractiveOutput]) -> None:

        message = args.message

        if message.find_cavities is not None:


            pass

        raise ValueError("Message does not contain a valid command.")

class Cavities(Module):

    def interactive_mode(self, options: Dict[str, Any]):
        pass

    def main(self, context: Context) -> Result:
        options = docopt(__doc__)
        
        if options.get('interactive'):
            return self.interactive_mode(options)