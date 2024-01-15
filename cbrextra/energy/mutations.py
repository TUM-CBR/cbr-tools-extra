from Bio.Data import IUPACData
import json
from os import path
import re
from typing import List, NamedTuple, TextIO

from . import sequences

ROTAMER_RE = re.compile(r"rot(?P<rot>\d+)")

POSITON_RE = re.compile(r"(resi\s+)?(?P<position>\d+)")

RESIDUES_3TO1 = dict(
    (key.upper(), value.upper())
    for key,value in IUPACData.protein_letters_3to1.items()
)

def three_to_one(three_letter_code : str) -> str:

    if len(three_letter_code) == 1:
        return three_letter_code.upper()

    return RESIDUES_3TO1.get(three_letter_code.upper(), "X")

class Mutation(NamedTuple):
    position : int
    mutation : str
    original : str

    @property
    def mutation_str(self) -> str:
        return "%s%i%s" % (self.original, self.position, self.mutation)

    @staticmethod
    def from_json_dict(
        json_dict : dict,
        sequence : str
    ) -> 'Mutation':

        mutation_str = json_dict['selection']
        position_match = POSITON_RE.match(mutation_str)

        if position_match is None:
            raise Exception("The string '%s' is not a valid mutation definition." % mutation_str)
        
        position = int(position_match.group('position'))

        return Mutation(
            position=position,
            mutation=three_to_one(json_dict['mutation']),
            original=sequence[position-1]
        )

class MutationResult(NamedTuple):
    structure_name : str
    mutations : List[Mutation]
    mutated_structure_name : str
    i_rotamer : int

    @property
    def mutation_str(self):
        return ",".join(m.mutation_str for m in self.mutations)

    @staticmethod
    def from_json_dict(json_dict : dict) -> 'MutationResult':

        name = json_dict['name']
        seq = sequences.get_sequence(name)
        mutations = [Mutation.from_json_dict(m, seq) for m in json_dict['mutations']]

        return MutationResult(
            structure_name=name[0:4],
            mutations=mutations,
            mutated_structure_name=name,
            i_rotamer=json_dict['rotamer']
        )

    @staticmethod
    def from_stream(stream : TextIO) -> 'MutationResult':
        return MutationResult.from_json_dict(json.load(stream))

    def get_energy_log_file(self, directory : str):
        return path.join(directory, "%s.log" % self.structure_name)

