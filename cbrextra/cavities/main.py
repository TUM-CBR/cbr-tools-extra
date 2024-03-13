"""Usage:
    cbrtools cavities interactive --input-points=<input_points_json> --radii-scale=<radii_scale> --empty-treshold=<empty-treshlod> --max_size_multiplier=<max_size_multiplier>
"""

from docopt import docopt
import json
import numpy as np
import sys
from typing import Any, Dict, TextIO

from cbrextra.cavities.data import InteractiveInput
from cbrextra.core.interactive.shared import ParseMessageArgs, SerializeMessageArgs

from ..core.interactive import *
from ..core.interactive.shared import OnMessageArgs
from ..core.module import Context, Module, Result
from . import algorithms
from .data import *

class CavitiesInstance:

    def __init__(
        self,
        points: Points,
        settings: Options
    ):
        self.__id = points.id
        self.__points = algorithms.to_spheres_cloud(points.points, points.radii, settings.radii_scale_or_default())
        self.__graph = algorithms.find_cavities(np.array(points.points), self.__points, settings)

    @property
    def id(self) -> str:
        return self.__id
    
    def get_cavities(self, min_volume: int = 2, max_volume: int = 125):
        cavities = self.__graph.get_cavities(min_volume, max_volume)

        return [cavity.to_cavity_model() for cavity in cavities]

class CavititesInteractive(InteractiveSpec[InteractiveInput, InteractiveOutput]):

    def __init__(
        self,
        points_set: List[Points],
        options: Options
    ):
        self.__cavities = [
            CavitiesInstance(points, options)
            for points in points_set
        ]

    def __handle_find_cavities(self, reqs: List[FindCavitiesArgs]):

        return CavitiesResult(
            cavities = {
                req.points_id: self.__handle_find_cavity(req)
                for req in reqs
            }
        )

    def __handle_find_cavity(self, req: FindCavitiesArgs):

        for cavity in self.__cavities:

            if cavity.id == req.points_id:
                return cavity.get_cavities(req.min_volume, req.max_volume)
            
        raise ValueError(f"The id '{req.points_id}' does not match any object.")
    
    def parse_message(self, args: ParseMessageArgs[InteractiveInput, InteractiveOutput]) -> InteractiveInput:
        return InteractiveInput(**args.payload)
    
    def serialize_message(self, args: SerializeMessageArgs[InteractiveInput, InteractiveOutput]) -> Dict[Any, Any]:
        return args.value.model_dump()

    def on_message(self, args: OnMessageArgs[InteractiveInput, InteractiveOutput]) -> None:

        message = args.message

        if message.find_cavities is not None:
            args.send(
                InteractiveOutput(
                    cavities_result=self.__handle_find_cavities(message.find_cavities)
                )
            )
            return 

        raise ValueError("Message does not contain a valid command.")

class Cavities(Module):

    def interactive_mode(
        self,
        options: Dict[str, Any],
        input_stream: TextIO = sys.stdin,
        output_stream: TextIO = sys.stdout
    ):

        with open(options['--input-points'], 'r') as points_json:
            points = [
                Points(**value)
                for value in json.load(points_json)
            ]

        run_interactive(
            CavititesInteractive(points, Options.from_args(options)),
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
    
module = Cavities()