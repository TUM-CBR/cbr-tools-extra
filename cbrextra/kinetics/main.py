"""Usage:
    cbrtools kinetics eval
    cbrtools kinetics interactive
"""

from docopt import docopt
import json
import sys
from typing import Any, Dict, TextIO

from ..core.module import Context, Module, Result
from ..core.interactive import *
from .compute import eval_model
from .data import *

class EnzymeKineticsInteractive(InteractiveSpec[InteractiveInput, InteractiveOutput]):

    def parse_message(self, args: 'ParseMessageArgs[InteractiveInput, InteractiveOutput]') -> InteractiveInput:
        return InteractiveInput(**args.payload)

    def serialize_message(self, args: 'SerializeMessageArgs[InteractiveInput, InteractiveOutput]') -> dict:
        return args.value.model_dump()

    def on_message(self, args: 'OnMessageArgs[InteractiveInput, InteractiveOutput]') -> None:

        message = args.message

        if message.eval_model is not None:
            result = eval_model(message.eval_model)
            args.send(
                InteractiveOutput(
                    eval_result = result
                )
            )
            return
        else:
            raise ValueError("Message contains no interactive command.")
        

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
        run_interactive(
            EnzymeKineticsInteractive(),
            input_stream,
            output_stream
        )
        return Result.success()

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if options.get('eval'):
            return self.eval_model(options)

        if options.get('interactive'):
            return self.interactive_mode(options)

        return Result.not_requested()

module = EnzymeKinetics()