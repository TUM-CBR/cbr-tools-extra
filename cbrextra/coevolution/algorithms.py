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

        msa_df = pd.merge(msa_df, msa_df, how='outer', suffixes=[LSUFFIX, RSUFFIX], 
        msa_df = msa_df[msa_df[K_POSITION_1] != msa_df[K_POSITION_2]]