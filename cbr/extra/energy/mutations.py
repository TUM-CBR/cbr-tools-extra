import json
from os import path
import re
from typing import List, NamedTuple, Self, TextIO

ROTAMER_RE = re.compile(r"rot(?P<rot>\d+)")

POSITON_RE = re.compile(r"(resi\s+)?(?P<position>\d+)")

class Mutation(NamedTuple):
    position : int
    mutation : str
    original : str

    @property
    def mutation_str(self : Self) -> str:
        return "%s%i%s" % (self.original, self.position, self.mutation)

    @staticmethod
    def from_json(
        json_dict : dict,
        sequence : str
    ) -> 'Mutation':

        mutation_str = json_dict['position']
        position_match = POSITON_RE.match(mutation_str)

        if position_match is None:
            raise Exception("The string '%s' is not a valid mutation definition." % mutation_str)
        
        position = int(position_match.group('position'))

        return Mutation(
            position=position,
            mutation=json_dict['mutation'],
            original=sequence[position]
        )

class MutationResult(NamedTuple):
    structure_name : str
    mutations : List[Mutation]
    original_residue : str
    mutation_position : int
    new_reside : str
    mutated_structure_name : str
    i_rotamer : int

    def to_json_dict(self):
        return {
            'structure_name': self.structure_name,
            'original_residue': self.original_residue,
            'mutation_position': self.mutation_position,
            'new_reside': self.new_reside,
            'mutated_structure_name': self.mutated_structure_name
        }



    @property
    def mutation_str(self):
        return ",".join(m.mutation_str for m in self.mutations)

    @staticmethod
    def from_json_dict(json_dict : dict) -> 'MutationResult':
        return MutationResult(**json_dict)

    @staticmethod
    def from_stream(stream : TextIO) -> 'MutationResult':
        return MutationResult.from_json_dict(json.load(stream))

    def get_energy_log_file(self, directory : str):
        return path.join(directory, "%s.log" % self.structure_name)

