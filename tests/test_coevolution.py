from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from typing import NamedTuple
from cbrextra import testutils

from cbrextra.coevolution.algorithms import CoEvolutionAnalysis

class CoevolutionTestSpec(NamedTuple):
    alignment_file: str

    def open_alignment(self) -> MultipleSeqAlignment:
        file = testutils.get_resource_path(
            __file__,
            "coevolution",
            self.alignment_file
        )
        return AlignIO.read(file, format='fasta')


SPEC_FORMYL_TRANSFERASE = CoevolutionTestSpec(
    alignment_file="FormylTransferase.aln.fasta"
)

class TestCoevolution:

    def test_open_alignment(self):

        spec = SPEC_FORMYL_TRANSFERASE

        CoEvolutionAnalysis.create(spec.open_alignment())

        assert True