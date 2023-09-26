from typing import Dict, Iterable, Optional, Tuple

from .data import DesignPrimersArgs, DesignPrimersResult, DesignPrimersResults, NamedTuple, PrimerOrganism, PrimerResult
from .melting_temp import MeltingTemp

CODON_SIZE = 3

DEFAULT_DNA_CONC = 500
DEFAULT_DNTP_CONC = 0.2

class Operations:
    class DesignPrimers(NamedTuple):
        operations : 'Operations'
        args : DesignPrimersArgs

        @property
        def organism(self) -> PrimerOrganism:
            return self.args.organism

        @property
        def start(self) -> int:
            return self.args.start

        @property
        def codon_count(self) -> int:
            return self.args.codon_count

        @property
        def min_length(self) -> int:
            return self.args.min_length

        @property
        def max_length(self) -> int:
            return self.args.max_length

        @property
        def sequence(self) -> str:
            return self.args.sequence

        @property
        def __step(self):
            return CODON_SIZE

        @property
        def codons(self) -> Dict[str, str]:
            return CODONS_MAP[self.organism]

        def design_primers(self) -> DesignPrimersResults:

            codon_count = self.start + self.codon_count*self.__step + self.__step

            for i in range(self.start, codon_count, self.__step):
                yield self.design_primer_at(i)

        def design_primer_at(self, position : int) -> DesignPrimersResult:
            tm_calc = self.operations.tm_calc
            results : DesignPrimersResult = {}

            for aa,codon in self.codons.items():

                def get_primers_for_codon() -> Iterable[PrimerResult]:
                    for (p_left, o_codon, p_right) in self.generate_primers_at(position):
                        tm_left = tm_calc.oligo_tm(p_left)
                        tm_right = tm_calc.oligo_tm(p_right)

                        seq = p_left + o_codon + p_right
                        tm_all = tm_calc.oligo_tm(seq)

                        yield PrimerResult(
                            left_primer=p_left,
                            tm_left=tm_left,
                            right_primer=p_right,
                            tm_right=tm_right,
                            inner_seq=codon,
                            tm_all=tm_all
                        )

                results[aa] = list(get_primers_for_codon())

            return results

        def __is_unique(self, primer_candidate: str):
            start_ix = 0

            while(start_ix < len(self.sequence)):

                next_ix = self.sequence.find(primer_candidate, start_ix)
                if next_ix >= 0 and start_ix > 0:
                    return False
                elif next_ix >= 0:
                    start_ix = next_ix + 1
                elif start_ix == 0:
                    raise Exception("The string %s does not occur in %s" % (primer_candidate, self.sequence))
                else:
                    return True

            return True

        def generate_primers_at(
            self,
            position : int,
        ) -> Iterable[Tuple[str, str, str]]:
            count = self.__step
            min_length = self.min_length
            max_length = self.max_length + 1

            for i in range(min_length, max_length):
                for j in range(min_length, max_length):
                    l_start = position - i
                    r_start = position + count
                    r_end = position + count + j

                    if l_start > 0 and r_start < len(self.sequence):

                        l_seq = self.sequence[l_start:position]
                        r_seq = self.sequence[r_start:r_end]
                        if self.__is_unique(l_seq) and self.__is_unique(r_seq):
                            yield (
                                l_seq,
                                self.sequence[position:r_start],
                                r_seq
                            )

    def __init__(self, tm_calc : Optional[MeltingTemp] = None):
        self.__tm_calc = \
            tm_calc or \
            MeltingTemp(
                dna_conc=DEFAULT_DNA_CONC,
                dntp_conc=DEFAULT_DNTP_CONC
            )

    @property
    def tm_calc(self) -> MeltingTemp:
        return self.__tm_calc

    def __design_primers(
        self,
        args: DesignPrimersArgs
    ) -> DesignPrimersResults:

        return Operations.DesignPrimers(
            args = args,
            operations=self
        ).design_primers()

    @staticmethod
    def design_primers(args: DesignPrimersArgs) -> DesignPrimersResults:
        return Operations().__design_primers(args)

P_PASTORIS_CODONS = {
    "Ala": "GCT",
    "Cys": "TGT",
    "Asp": "GAC",
    "Glu": "GAG",
    "Phe": "TTC",
    "Gly": "GGT",
    "His": "CAC",
    "Ile": "ATT",
    "Lys": "AAG",
    "Leu": "TTG",
    "Met": "ATG",
    "Asn": "AAC",
    "Pro": "CCA",
    "Gln": "CAA",
    "Arg": "AGA",
    "Ser": "TCT",
    "Thr": "ACT",
    "Val": "GTT",
    "Trp": "TGG",
    "Tyr": "TAC"
}

E_COLI_CODONS = {
    "Ala": "GCG",
    "Cys": "TGC",
    "Asp": "GAC",
    "Glu": "GAA",
    "Phe": "TTC",
    "Gly": "GGT",
    "His": "CAC",
    "Ile": "ATC",
    "Lys": "AAA",
    "Leu": "CTG",
    "Met": "ATG",
    "Asn": "AAC",
    "Pro": "CCG",
    "Gln": "CAG",
    "Arg": "CGT",
    "Ser": "TCT",
    "Thr": "ACC",
    "Val": "GTT",
    "Trp": "TGG",
    "Tyr": "TAC"
}

CODONS_MAP : Dict[PrimerOrganism, Dict[str, str]] = {
    'E_COLI': E_COLI_CODONS,
    'P_PASTORIS': P_PASTORIS_CODONS
}