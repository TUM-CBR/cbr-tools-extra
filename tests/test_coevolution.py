from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from typing import Any, NamedTuple, cast
from cbrextra import testutils

from cbrextra.coevolution.data import Scoring
from cbrextra.coevolution.algorithms import CoEvolutionAnalysis, K_POSITION, K_RESIDUE, K_SEQUENCE, RESIDUE_TO_INT

class CoevolutionTestSpec(NamedTuple):
    alignment_file: str

    def open_alignment(self) -> MultipleSeqAlignment:
        file = testutils.get_resource_path(
            __file__,
            "coevolution",
            self.alignment_file
        )
        return AlignIO.read(file, format='fasta')


TEST_SPECS = [
    CoevolutionTestSpec(
        alignment_file="Test1.aln.fasta"
    ),
    CoevolutionTestSpec(
        alignment_file="Test2.aln.fasta"
    )
]

class TestCoevolution:

    def test_open_alignment(self):

        for spec in TEST_SPECS:

            msa = spec.open_alignment()
            analysis = CoEvolutionAnalysis.create(msa)

            pos = 123
            seq = 123

            data = analysis.data
            record = data[(data[K_POSITION] == pos) & (data[K_SEQUENCE] == seq)]
            expected_resi = RESIDUE_TO_INT[msa[seq]._seq[pos].upper()]
            actual_resi = list(cast(Any, record[K_RESIDUE]))[0]

            assert expected_resi == actual_resi

    def test_score_alignment(self):

        for spec in TEST_SPECS:
            
            msa = spec.open_alignment()
            analysis = CoEvolutionAnalysis.create(msa)

            position = 71
            results_per_position = 20


            scores = analysis.query_scores(
                [position],
                results_per_position=results_per_position,
                scoring=Scoring()
            )

            assert True