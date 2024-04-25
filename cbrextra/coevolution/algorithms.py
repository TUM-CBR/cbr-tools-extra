from Bio.Align import MultipleSeqAlignment
from Bio.Data.IUPACData import protein_letters
import pandas as pd
from typing import cast, Iterable, List, NamedTuple, Optional

from .data import *

RESIDUE_TO_INT = {
    res: i
    for i,res in enumerate(protein_letters)
}

IGAP = -1
RESIDUE_TO_INT['-'] = IGAP

LSUFFIX = '_res1'
RSUFFIX = '_res2'
K_SEQUENCE = 'sequence'
K_POSITION = 'position'
K_RESIDUE = 'residue'
K_POSITION_1 = K_POSITION + LSUFFIX
K_POSITION_2 = K_POSITION + RSUFFIX
K_RESIDUE_1 = K_RESIDUE + LSUFFIX
K_RESIDUE_2 = K_RESIDUE + RSUFFIX
K_OCCURRENCE_COUNT = 'occurence count'
K_OCCURENCE_SINGLE_COUNT = 'occurence single count'
K_HETEROGENEITY_SCORE = "occurence heterogeneity"
K_SCORE_ALL = "score"

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

        for seq_id, seq_record in enumerate(alignment):
            for (pos, res) in enumerate(seq_record._seq):
                sequence[i] = seq_id
                position[i] =  pos
                residue[i] = RESIDUE_TO_INT[res.upper()]
                i += 1

        msa_df = pd.DataFrame(
            data = {
                K_SEQUENCE: sequence,
                K_POSITION: position,
                K_RESIDUE: residue
            }
        )

        return CoEvolutionAnalysis(
            alignment = alignment,
            data = msa_df
        )
    
    @classmethod
    def score_conserved(cls, result: pd.DataFrame) -> 'pd.Series[float]':
        """
        This scoring amis to exclude positions where residues are simply
        always conserved. These are potenitally important relations but
        not for the purpose of co-evolution. The idea is that the closer
        the residue pair is to being fully conserved accross the whole
        set of sequences, the more the pair gets penalized.

        Parameters
        ----------
        result: pd.DataFrame
            A DataFrame containing the K_OCCURRENCE_SCORE column with the occurence
            or residue pairs relative to the total.
        """

        series = result[K_OCCURRENCE_SCORE]
        low_mask = result[K_OCCURRENCE_SCORE] > 0.5
        penalty = float(2)*(series - float(0.5)) * low_mask
        return float(1) - penalty
        
    def score_occurence(self, result: pd.DataFrame) ->  'pd.Series[float]':
        """
        This function scores the occurence of a specific pair of residues relative
        to the total number of sequences in the alignment.

        
        Parameters
        ----------
        series: pd.DataFrame
            A DataFrame containing the counts per residue pair as the column named
            K_OCCURRENCE_COUNT.
        """

        return result[K_OCCURRENCE_COUNT] / len(self.alignment)
    
    def score_symmetry(self, result: pd.DataFrame) -> 'pd.Series[float]':

        res_1 = result[K_RESIDUE_1].to_numpy()
        pos_1 = result[K_POSITION_1].to_numpy()

        res_2 = result[K_RESIDUE_2].to_numpy()
        pos_2  = result[K_POSITION_2].to_numpy()

        occurrence = result[K_OCCURRENCE_SCORE]

        K_LEFT = "occurrence left"
        K_RIGHT = "occurrence right"

        compare_df = pd.merge(
            pd.DataFrame(
                data={
                    K_RESIDUE_1: res_1,
                    K_POSITION_1: pos_1,
                    K_RESIDUE_2: res_2,
                    K_POSITION_2: pos_2,
                    K_LEFT: occurrence
                }
            ),
            pd.DataFrame(
                data={
                    K_RESIDUE_1: res_2,
                    K_POSITION_1: pos_1,
                    K_RESIDUE_2: res_1,
                    K_POSITION_2: pos_2,
                    K_RIGHT: occurrence
                }
            ),
            how='left',
            on=[K_RESIDUE_1, K_POSITION_1, K_RESIDUE_2, K_POSITION_2]
        ).fillna({
            K_LEFT: 0,
            K_RIGHT: 0
        })

        left = compare_df[K_LEFT]
        right = compare_df[K_RIGHT]
        score = float(1) - ((left - right).abs() / (left + right))

        return score

    @classmethod
    def score_exclusivity(cls, result: pd.DataFrame) -> 'pd.Series[float]':
        """
        This scoring aims to capture how often does this specifi residue
        pairs occurs relative to the first residue of the pair occuring
        paired with other residues. The intuition is that if two residues
        are co-evolving, then one should expect to observe the second residue
        whenever the first residue occurs. A score of 1 means that every
        time the first residue of the pair occurs at the specified position,
        then the second residue of the pair will follow.

        Parameters
        ----------
        result: pd.DataFrame
            A pandas DataFrame containing the colum K_OCCURRENCE_COUNT with
            the counts of occurence of the reside pair and a column called
            K_OCCURENCE_SINGLE_COUNT containing the counts of the occurence
            of the first residue of the pair paired with any other residue.
        """
        return result[K_OCCURRENCE_COUNT] / result[K_OCCURENCE_SINGLE_COUNT]
    
    def __get_paired_dataframe(self, positions: List[int]) -> pd.DataFrame:
        """
        Construct the DataFrame that contains the positions of interest paired
        with every other position in the alignment (except itself).
        """
        if len(positions) < 1:
            raise ValueError(f"{positions}: The length must be at least 1.")
        

        position = positions[0]
        data = self.data
        mask = data[K_POSITION] == position

        for position in positions[1:]:
            mask |= data[K_POSITION] == position

        data = pd.merge(data[mask], data, how='outer', suffixes=[LSUFFIX, RSUFFIX], on=[K_SEQUENCE])

        return data[data[K_POSITION_1] != data[K_POSITION_2]]

    def score_positions(
        self,
        positions: List[int],
        scoring: Scoring
    ) -> pd.DataFrame:
        """
        Query the co-evolution analysis using the correlated occurence of the
        residues at the positions provided with residues in any other positions
        of the sequence alignment.
        """

        data =  self.__get_paired_dataframe(positions)
        mask_gaps = (data[K_RESIDUE_1] != IGAP) & (data[K_RESIDUE_2] != IGAP)
        mask_same = (data[K_RESIDUE_1] != data[K_RESIDUE_2])
        mask_for_pairs = mask_gaps & mask_same
        pair_counts = data[mask_for_pairs].groupby(by=[K_POSITION_1, K_POSITION_2, K_RESIDUE_1, K_RESIDUE_2]).count()
        pair_counts_index = pair_counts.index.to_frame(index=False)
        single_counts = data.groupby(by=[K_POSITION_1, K_POSITION_2, K_RESIDUE_1]).count()
        single_counts_index = single_counts.index.to_frame(index=False)

        result = pd.DataFrame(
            data = {
                K_POSITION_1: pair_counts_index[K_POSITION_1].to_numpy(),
                K_POSITION_2: pair_counts_index[K_POSITION_2].to_numpy(),
                K_RESIDUE_1: pair_counts_index[K_RESIDUE_1].to_numpy(),
                K_RESIDUE_2: pair_counts_index[K_RESIDUE_2].to_numpy(),
                K_OCCURRENCE_COUNT: pair_counts[K_SEQUENCE].to_numpy()
            }
        )

        result = pd.merge(
            result,
            pd.DataFrame(
                data = {
                    K_POSITION_1: single_counts_index[K_POSITION_1].to_numpy(),
                    K_POSITION_2: single_counts_index[K_POSITION_2].to_numpy(),
                    K_RESIDUE_1: single_counts_index[K_RESIDUE_1].to_numpy(),
                    K_OCCURENCE_SINGLE_COUNT: single_counts[K_SEQUENCE].to_numpy()
                }
            ),
            on=[K_POSITION_1, K_POSITION_2, K_RESIDUE_1],
            how='left'
        )

        result[K_OCCURRENCE_SCORE] = score_occurence = self.score_occurence(result)
        result[K_EXLUSIVITY_SCORE] = score_exclusivity = result[K_OCCURRENCE_SCORE] / result[K_OCCURENCE_SINGLE_COUNT]
        result[K_CONSERVED_SCORE] = score_conserved = self.score_conserved(result)
        result[K_SYMMETRY_SCORE] = score_symmetry = self.score_symmetry(result)
        result[K_SCORE_ALL] = (
            scoring.occurence_weight * score_occurence
            + scoring.exclusivity_weight * score_exclusivity
            + scoring.conserved_weight * score_conserved
            + scoring.symmetry_weight * score_symmetry
        ) / (scoring.occurence_weight + scoring.exclusivity_weight + scoring.conserved_weight + scoring.symmetry_weight)

        return result
    
    @classmethod
    def __to_position_result(cls, ranked_position_scores: pd.DataFrame) -> CoevolutionPosition:
        positions: 'pd.Series[int]' = ranked_position_scores[K_POSITION_1] 
        position = positions.iat[0]

        if not (positions == position).all():
            raise ValueError()


        positions_2: 'pd.Series[int]' = ranked_position_scores[K_POSITION_2]
        residues_1: 'pd.Series[int]' = ranked_position_scores[K_RESIDUE_1]
        residues_2: 'pd.Series[int]' = ranked_position_scores[K_RESIDUE_2]
        scores: 'pd.Series[float]' = ranked_position_scores[K_SCORE_ALL]
        score_occurence: 'pd.Series[float]' = ranked_position_scores[K_OCCURRENCE_SCORE]
        score_exclusivity: 'pd.Series[float]' = ranked_position_scores[K_EXLUSIVITY_SCORE]
        score_conserved: 'pd.Series[float]' = ranked_position_scores[K_CONSERVED_SCORE]
        score_symmetry: 'pd.Series[float]' = ranked_position_scores[K_SYMMETRY_SCORE]

        return CoevolutionPosition(
            position = position,
            by_position = {
                positions_2[i]: CoevolutionEntry(
                    residue_1=protein_letters[residues_1[i]],
                    residue_2=protein_letters[residues_2[i]],
                    score=scores[i],
                    score_occurence=score_occurence[i],
                    score_exclusivity=score_exclusivity[i],
                    score_conserved=score_conserved[i],
                    score_symmetry=score_symmetry[i]
                )
                for i in cast(Iterable[int], ranked_position_scores.index)
            }
        )

    def query_scores(
        self,
        positions: List[int],
        results_per_position: int,
        scoring: Scoring
    ) -> CoevolutionResults:
        scores = self.score_positions(positions, scoring)

        coevolution_positions = {
            position : self.__to_position_result(score)
            for position in positions
            for score in [scores[scores[K_POSITION_1] == position]]
            for ranks in [score[K_SCORE_ALL].rank(ascending=False)]
            for score in [score[ranks <= results_per_position]]
        }

        return CoevolutionResults(
            positions = coevolution_positions
        )