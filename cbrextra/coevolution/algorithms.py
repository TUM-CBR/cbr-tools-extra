from Bio.Align import MultipleSeqAlignment
from Bio.Data.IUPACData import protein_letters
import pandas as pd
from typing import List, NamedTuple, Optional

RESIDUE_TO_INT = {
    res: i
    for i,res in enumerate(protein_letters)
}

LSUFFIX = '_seq1'
RSUFFIX = '_seq2'
K_SEQUENCE = 'sequence'
K_POSITION = 'position'
K_RESIDUE = 'residue'
K_POSITION_1 = K_POSITION + LSUFFIX
K_POSITION_2 = K_POSITION + RSUFFIX
K_RESIDUE_1 = K_RESIDUE + LSUFFIX
K_RESIDUE_2 = K_RESIDUE + RSUFFIX
K_OCCURRENCE_COUNT = 'occurence count'
K_OCCURENCE_SINGLE_COUNT = 'occurence single count'
K_OCCURRENCE_SCORE = 'occurrence'
K_EXLUSIVITY_SCORE = 'exclusivity'

class CoEvolutionAnalysis(NamedTuple):
    """
    Class responsible for storing the data contained in a
    mustiple sequence analysis in such a way it is efficient for
    performing co-evolution analysis.

    Attributes
    ----------
    alignment: MultipleSeqAlignment
        The alignment from which this class was derived
    data: pd.DataFrame
        The data of the alignment necessary for co-evolution analysis
        stored as a pandas data frame for efficient retrival.
    """

    alignment: MultipleSeqAlignment
    data: pd.DataFrame

    @classmethod
    def create(cls, alignment: MultipleSeqAlignment) -> 'CoEvolutionAnalysis':
        """
        Create an instance  of this  class by correctly indexing the data of the
        analysis. This consists of fist storing the sequence id, residue  position and
        residue name in a data frame. We then compute the cross product of this DataFrame
        with itself to have all the possible pair of residues
        """
        size = alignment.get_alignment_length() * len(alignment)
        sequence: List[Optional[int]] = [None for _ in range(size)]
        position: List[Optional[int]] = [None for _ in range(size)]
        residue: List[Optional[int]] = [None for _ in range(size)]
        i = 0

        for seq_id, sequence in enumerate(alignment):
            for (pos, res) in enumerate(sequence._seq):
                sequence[i] = seq_id
                position[i] =  pos
                residue[i] = RESIDUE_TO_INT[res.upper()]

        msa_df = pd.DataFrame(
            data = {
                K_SEQUENCE: sequence,
                K_POSITION: position,
                K_RESIDUE: residue
            }
        )

        msa_df = pd.merge(msa_df, msa_df, how='outer', suffixes=[LSUFFIX, RSUFFIX], on=[K_SEQUENCE])

        return CoEvolutionAnalysis(
            alignment = alignment,
            data = msa_df[msa_df[K_POSITION_1] != msa_df[K_POSITION_2]]
        )
    
    def score_positions(
        self,
        positions: List[int],
        center: float
    ):
        """
        Query the co-evolution analysis using the correlated occurence of the
        residues at the positions provided with residues in any other positions
        of the sequence alignment
        """

        assert len(positions) > 0, "Position has to be provided to query."

        position = positions[0]
        data = self.data
        mask = data[K_POSITION_1] == position

        for position in positions[1:]:
            mask |= data[K_POSITION_1]  == position

        data = data[mask]
        pair_counts = data.groupby(by=[K_POSITION_1, K_POSITION_2, K_RESIDUE_1, K_RESIDUE_2]).count()
        pair_counts_index = pair_counts.index.to_frame(index=False)
        single_counts = data.groupby(by=[K_POSITION_1, K_POSITION_2, K_RESIDUE_1]).count()
        single_counts_index = single_counts.index.to_frame(index=False)

        result = pd.DataFrame(
            data = {
                K_POSITION_1: pair_counts_index[K_POSITION_1],
                K_POSITION_2: pair_counts_index[K_POSITION_2],
                K_RESIDUE_1: pair_counts_index[K_RESIDUE_1],
                K_RESIDUE_2: pair_counts_index[K_RESIDUE_2],
                K_OCCURRENCE_COUNT: pair_counts[K_SEQUENCE]
            }
        )

        result = pd.merge(
            result,
            pd.DataFrame(
                data = {
                    K_POSITION_1: single_counts_index[K_POSITION_1],
                    K_POSITION_2: single_counts_index[K_POSITION_2],
                    K_RESIDUE_1: single_counts_index[K_RESIDUE_1],
                    K_OCCURENCE_SINGLE_COUNT: single_counts[K_OCCURENCE_SINGLE_COUNT]
                }
            ),
            on=[K_POSITION_1, K_POSITION_2, K_RESIDUE_1],
            how='left'
        )

        result[K_OCCURRENCE_SCORE] = result[K_OCCURRENCE_COUNT] / len(self.alignment)
        result[K_EXLUSIVITY_SCORE] = result[K_OCCURRENCE_SCORE] / result[K_OCCURENCE_SINGLE_COUNT]

