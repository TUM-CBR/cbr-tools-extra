"""Usage:
    cbrtools coevolution interactive --input-msa=<input-msa>
"""

from typing import Dict
from docopt import docopt
import json

from ..core.interactive import *
from ..core.interactive.shared import OnMessageArgs, ParseMessageArgs, SerializeMessageArgs

from .data import InteractiveInput, InteractiveOutput

class CoevolutionInteractive(InteractiveSpec[InteractiveInput, InteractiveOutput]):

    def parse_message(self, args: ParseMessageArgs[InteractiveInput, InteractiveOutput]) -> InteractiveInput:
        return InteractiveInput(**args.payload)
    
    def serialize_message(self, args: SerializeMessageArgs[InteractiveInput, InteractiveOutput]) -> Dict[Any, Any]:
        return args.value.model_dump()
    
    def on_message(self, args: OnMessageArgs[InteractiveInput, InteractiveOutput]) -> None:
        return super().on_message(args)