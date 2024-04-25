"""Usage:
    cbrtools coevolution interactive --input-msa=<input-msa>
"""

from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from docopt import docopt
from typing import Dict
import sys
from typing import cast, TextIO

from ..core.interactive import *
from ..core.interactive.shared import OnMessageArgs, ParseMessageArgs, SerializeMessageArgs
from ..core.module import Context, Module, Result

from .algorithms import CoEvolutionAnalysis
from .data import CoevolutionResults, InteractiveRequest, InteractiveResponse, Query

class CoevolutionInteractive(InteractiveSpec[InteractiveRequest, InteractiveResponse]):

    def __init__(
        self,
        alignment: MultipleSeqAlignment
    ) -> None:
        super().__init__()

        self.__coevolution = CoEvolutionAnalysis.create(alignment)

    def parse_message(self, args: ParseMessageArgs[InteractiveRequest, InteractiveResponse]) -> InteractiveRequest:
        return InteractiveRequest(**args.payload)
    
    def serialize_message(self, args: SerializeMessageArgs[InteractiveRequest, InteractiveResponse]) -> Dict[Any, Any]:
        return args.value.model_dump()
    
    def __handle_query(self, query: Query) -> CoevolutionResults:
        return self.__coevolution.query_scores(
            query.positions,
            query.max_results,
            query.scoring
        )

    def on_message(self, args: OnMessageArgs[InteractiveRequest, InteractiveResponse]) -> None:

        message = args.message

        if message.query is not None:
            args.send(
                InteractiveResponse(coevolution=self.__handle_query(message.query))
            )
        else:
            raise ValueError("Message contains no requests.")
        
class Coevolution(Module):

    def interactive_mode(
        self,
        options: Dict[str, Any],
        input_stream: TextIO = sys.stdin,
        output_stream: TextIO = sys.stdout
    ):
        
        msa = AlignIO.read(options['--input-points'], format='fasta')

        run_interactive(
            CoevolutionInteractive(cast(MultipleSeqAlignment, msa)),
            input_stream,
            output_stream
        )

    def main(
        self,
        context: Context,
        argv: Optional[List[str]] = None,
        input_stream: TextIO = sys.stdin,
        output_stream: TextIO = sys.stdout
    ) -> Result:
        options = docopt(__doc__, argv=argv)

        if options.get('interactive'):
            self.interactive_mode(options, input_stream=input_stream, output_stream=output_stream)
            return Result.success()
        
        return Result.not_requested()
    
module = Coevolution()