from Bio.Align import MultipleSeqAlignment
from Bio.Data.IUPACData import protein_letters
import numpy as np
import pandas as pd
from typing import cast, Iterable, List, NamedTuple, Optional

from .data import *

RESIDUE_TO_INT = {
    res: i
    for i,res in enumerate(protein_letters)
}
TOTAL_RESIDUES = len(protein_letters)

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
K_OCCURENCE_SINGLE_COUNT_1 = K_OCCURENCE_SINGLE_COUNT + LSUFFIX
K_OCCURENCE_SINGLE_COUNT_2 = K_OCCURENCE_SINGLE_COUNT + RSUFFIX
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
    
    def score_confidence(
        self,
        result: pd.DataFrame,
        treshold: float = 0.025,
        center: float = 0.5,
        concaveness: float = 8
    ) -> 'pd.Series[float]':
        """
        This function scores how much confidence one can have on the results of
        the analysis. The idea is that positions with too much or too little
        overall conservation are less reliable than positions with a modest
        degree of conservation.


        Parameters
        ----------
        result: pd.DataFrame
            A DataFrame containing the K_OCCURRENCE_SCORE column with the occurence
            or residue pairs relative to the total.
        treshold: The minimum occurrence a residue must have at a position in the
            MSA to be counted as present.
        center: When weighting the scores, the number of present residues is divided
            by the total number of residues (20) at each position. The result is then
            weighted with the function e^(-concaveness*(x - center)^2) to assign
            the confidence socre at each position in the MSA. This parameter controls
            the "center" of the funciton.
        concaveness: This controls the concaveness of the function mentioned above.
        """

        freqs = np.vstack([
            freq
            for k,freq in self.alignment.alignment.frequencies.items()
                if k != "-"
        ]).T / len(self.alignment)

        counts_normalized = np.count_nonzero(freqs > treshold, axis=1) / TOTAL_RESIDUES
        counts_scoped = np.exp(-concaveness*np.power(counts_normalized - center, 2))
        counts_df = pd.DataFrame(
                data={
                    K_CONFIDENCE_SCORE: counts_scoped,
                    K_POSITION: np.arange(0, self.alignment.get_alignment_length())
                })
        
        score_1 = K_CONFIDENCE_SCORE + "_1"
        score_2 = K_CONFIDENCE_SCORE + "_2"

        scores_df = pd.merge(
            result,
            counts_df.rename({K_POSITION: K_POSITION_1, K_CONFIDENCE_SCORE: score_1}, axis=1),
            on=[K_POSITION_1],
            how='left'
        )

        scores_df = pd.merge(
            scores_df,
            counts_df.rename({K_POSITION: K_POSITION_2, K_CONFIDENCE_SCORE: score_2}, axis=1),
            on=[K_POSITION_2],
            how='left'
        )

        return (scores_df[score_1] + scores_df[score_2]) / 2
        
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
    
    @classmethod
    def score_symmetry(cls, result: pd.DataFrame) -> 'pd.Series[float]':
        """
        Symmetry scoring identifies residues at two positions which
        appear in an alternating fashion. The idea is that co-evolving
        residues will appear together, regardless of whether they appear
        on one side or the other.
        """

        K_LEFT = "occurrence left"
        K_RIGHT = "occurrence right"
        
        compare_df = pd.merge(
            result.rename({K_OCCURRENCE_SCORE: K_LEFT}, axis=1),
            result.rename(
                {
                    K_RESIDUE_1: K_RESIDUE_2,
                    K_RESIDUE_2: K_RESIDUE_1,
                    K_OCCURRENCE_SCORE: K_RIGHT
                },
                axis=1
            ),
            on=[K_RESIDUE_1, K_POSITION_1, K_RESIDUE_2, K_POSITION_2],
            how='left'
        ).fillna({
            K_LEFT: 0,
            K_RIGHT: 0
        })

        # If the same resiude appears in the two positions,
        # it always has a score of 1. We set it to 0.5 as
        # there is no way to know if there is symmetry or not
        identical_mask = (result[K_RESIDUE_1] == result[K_RESIDUE_2]) * 0.5

        left = compare_df[K_LEFT]
        right = compare_df[K_RIGHT]
        score = float(1) - ((left - right).abs() / (left + right))


        return score - identical_mask

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
        return 2 * result[K_OCCURRENCE_COUNT] / (result[K_OCCURENCE_SINGLE_COUNT_1] + result[K_OCCURENCE_SINGLE_COUNT_2])
    
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

        data = self.__get_paired_dataframe(positions)
        mask_gaps = (data[K_RESIDUE_1] != IGAP) & (data[K_RESIDUE_2] != IGAP)
        #mask_same = (data[K_RESIDUE_1] != data[K_RESIDUE_2])
        mask_for_pairs = mask_gaps # & mask_same
        pair_counts = data[mask_for_pairs].groupby(by=[K_POSITION_1, K_POSITION_2, K_RESIDUE_1, K_RESIDUE_2]).count()
        pair_counts_index = pair_counts.index.to_frame(index=False)

        singles = self.data.groupby(by=[K_POSITION, K_RESIDUE]).count()
        singles_index = singles.index.to_frame(index=False)
        singles = pd.DataFrame(
            data = {
                K_POSITION: singles_index[K_POSITION].to_numpy(),
                K_RESIDUE: singles_index[K_RESIDUE].to_numpy(),
                K_OCCURRENCE_COUNT: singles[K_SEQUENCE].to_numpy()
            }
        )

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
            singles.rename(
                {
                    K_POSITION: K_POSITION_1,
                    K_RESIDUE: K_RESIDUE_1,
                    K_OCCURRENCE_COUNT: K_OCCURENCE_SINGLE_COUNT_1
                },
                axis=1
            ),
            on=[K_POSITION_1, K_RESIDUE_1],
            how='left'
        )

        result = pd.merge(
            result,
            singles.rename(
                {
                    K_POSITION: K_POSITION_2,
                    K_RESIDUE: K_RESIDUE_2,
                    K_OCCURRENCE_COUNT: K_OCCURENCE_SINGLE_COUNT_2
                },
                axis=1
            ),
            on=[K_POSITION_2, K_RESIDUE_2],
            how='left'
        )

        result[K_OCCURRENCE_SCORE] = score_occurrence = self.score_occurence(result)
        result[K_EXLUSIVITY_SCORE] = score_exclusivity = self.score_exclusivity(result)
        result[K_CONFIDENCE_SCORE] = score_confidence = self.score_confidence(
            result,
            treshold=scoring.confidence_treshold,
            center=scoring.confidence_center,
            concaveness=scoring.confidence_concaveness
        )
        result[K_SYMMETRY_SCORE] = score_symmetry = self.score_symmetry(result)
        result[K_SCORE_ALL] = (
            + scoring.exclusivity_weight_scaled * score_exclusivity
            + scoring.confidence_weight_scaled * score_confidence
            + scoring.symmetry_weight_scaled * score_symmetry
            + scoring.occurrence_weight_scaled * score_occurrence
        )

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
        score_confidence: 'pd.Series[float]' = ranked_position_scores[K_CONFIDENCE_SCORE]
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
                    score_confidence=score_confidence[i],
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
            for ranks in [score[K_SCORE_ALL].rank(ascending=False, method='first')]
            for score in [score[ranks <= results_per_position]]
        }

        return CoevolutionResults(
            positions = coevolution_positions
        )