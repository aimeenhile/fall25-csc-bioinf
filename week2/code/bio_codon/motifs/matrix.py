# bio_codon/motifs/matrix.py

"""Support for various forms of sequence motif matrices.

Implementation of frequency (count) matrices, position-weight matrices,
and position-specific scoring matrices.
"""

import math
from typing import TypeVar, Generic, list, dict, Optional, Any

# --- Constants & Helpers ---
DNA_ALPHABET: list[str] = ['A', 'C', 'G', 'T']
TMatrix = TypeVar('TMatrix', float, int)

def _safe_log2(x: float, background: float) -> float:
    """Calculates log2(x / background), handling log(0)."""
    if x <= 0.0:
        return -math.inf
    if background <= 0.0:
        # Should not happen with valid background
        return math.inf 
        
    return math.log2(x) - math.log2(background)


# --- Base Matrix Class ---

class GenericPositionMatrix(Generic[TMatrix]):
    """Base class for position-specific matrices (PSM)."""
    matrix: list[dict[str, TMatrix]]
    length: int
    alphabet: list[str]

    def __init__(self, alphabet: list[str], values: list[dict[str, TMatrix]]):
        """Initializes the matrix with values."""
        self.alphabet = alphabet
        self.matrix = values
        self.length = len(values)
        
    def __len__(self) -> int:
        return self.length

    def __str__(self) -> str:
        """String representation of the matrix."""
        s: str = f"<{self.__class__.__name__} of length {self.length}>\n"
        header: str = "Pos\t" + "\t".join(self.alphabet)
        s += header + "\n"

        for i in range(self.length):
            pos_data = self.matrix[i]
            scores: list[str] = []
            for base in self.alphabet:
                val: TMatrix = pos_data.get(base, 0)
                # Format based on type (float for PWM/Freq, int for Counts)
                if isinstance(val, float):
                    scores.append(f"{val:.4f}")
                else:
                    scores.append(str(val))
            s += f"{i}\t{'\t'.join(scores)}\n"
        return s
        
    def __getitem__(self, index: int) -> dict[str, TMatrix]:
        """Allows access to a specific position (column)."""
        return self.matrix[index]

# --- Count Matrix ---

class CountMatrix(GenericPositionMatrix[int]):
    """Matrix containing raw counts of nucleotides."""
    
    def __init__(self, sequences: list[str]):
        """Initializes the CountMatrix from a list of aligned sequences."""
        if not sequences:
            super().__init__(DNA_ALPHABET, [])
            return

        L: int = len(sequences[0])
        # Initialize the matrix structure
        raw_counts: list[dict[str, int]] = [
            {base: 0 for base in DNA_ALPHABET} for _ in range(L)
        ]
        
        # Populate counts
        for seq in sequences:
            if len(seq) != L:
                print(f"Warning: Skipping unaligned sequence of length {len(seq)}")
                continue
            for i, base_char in enumerate(seq):
                base = str(base_char).upper()
                if base in DNA_ALPHABET:
                    raw_counts[i][base] += 1
        
        super().__init__(DNA_ALPHABET, raw_counts)

# --- Frequency Matrix ---

class FrequencyPositionMatrix(GenericPositionMatrix[float]):
    """Matrix where values are normalized frequencies."""
    
    def __init__(self, 
                 counts: CountMatrix, 
                 pseudocount: float = 0.5):
        """Calculates frequencies from a CountMatrix, applying pseudocounts."""
        
        freq_matrix_data: list[dict[str, float]] = []
        num_bases: int = len(DNA_ALPHABET)
        
        for pos_counts in counts.matrix:
            
            total_instances: int = 0
            for base in DNA_ALPHABET:
                total_instances += pos_counts.get(base, 0)
            
            total_after_pc: float = float(total_instances) + (pseudocount * num_bases)

            freq_pos: dict[str, float] = {}
            for base in DNA_ALPHABET:
                count: float = float(pos_counts.get(base, 0))
                freq: float = (count + pseudocount) / total_after_pc
                freq_pos[base] = freq
            
            freq_matrix_data.append(freq_pos)

        super().__init__(DNA_ALPHABET, freq_matrix_data)
        
    def consensus(self) -> str:
        """Returns the consensus sequence."""
        result: list[str] = []
        for pos_data in self.matrix:
            max_val: float = -1.0
            max_base: str = 'N'
            for base in DNA_ALPHABET:
                val: float = pos_data.get(base, 0.0)
                if val > max_val:
                    max_val = val
                    max_base = base
            result.append(max_base)
        return "".join(result)

# --- Position Weight Matrix ---

class PositionWeightMatrix(GenericPositionMatrix[float]):
    """Matrix of log-odds scores (PWM/PSSM)."""
    min_score: float
    max_score: float

    def __init__(self, 
                 freq_matrix: FrequencyPositionMatrix,
                 background: dict[str, float] = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}):
        """Calculates the PWM from a FrequencyPositionMatrix and background model."""
        
        pwm_matrix: list[dict[str, float]] = []
        
        for pos_data in freq_matrix.matrix:
            pwm_pos: dict[str, float] = {}
            for base in DNA_ALPHABET:
                freq: float = pos_data.get(base, 0.0)
                bg: float = background.get(base, 0.25)
                pwm_pos[base] = _safe_log2(freq, bg)
            pwm_matrix.append(pwm_pos)
            
        super().__init__(DNA_ALPHABET, pwm_matrix)
        
        # Calculate min/max scores
        min_s: float = 0.0
        max_s: float = 0.0
        for pos in self.matrix:
            if len(pos.values()) > 0:
                min_s += min(pos.values())
                max_s += max(pos.values())

        self.min_score = min_s
        self.max_score = max_s
        
    def calculate_score(self, sequence: str) -> float:
        """Calculates the score of a given sequence against the PWM."""
        if len(sequence) != self.length:
            return -math.inf

        score: float = 0.0
        
        for i, base_char in enumerate(sequence):
            base = str(base_char).upper()
            pos_data = self.matrix[i]
            
            # Use the score for the base, or the minimum score for an unknown base (N)
            position_min_score: float = min(pos_data.values()) if len(pos_data) > 0 else 0.0
            score += pos_data.get(base, position_min_score) 

        return score


class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """Class for the support of Position Specific Scoring Matrix calculations."""
    pass

    """
    def __init__(self, alphabet: str, pssm_matrix: Dict[str, List[float]]):
        """Initialize the PSSM."""
        super().__init__(alphabet, pssm_matrix)
        self.pssm = self.pwm # PSSM uses the same underlying structure

    @property
    def max(self) -> float:
        """Maximum possible score for the motif."""
        return self._array.max(axis=0).sum()

    @property
    def min(self) -> float:
        """Minimum possible score for the motif."""
        return self._array.min(axis=0).sum()

    def calculate(self, sequence: Any) -> np.ndarray:
        """Calculate the PSSM score for all positions in a sequence."""
        if not self.alphabet.is_superset(sequence.alphabet):
             raise ValueError("Sequence alphabet is not compatible with PSSM alphabet.")

        sequence_str = str(sequence)
        seq_len = len(sequence_str)
        motif_len = self.length
        
        if seq_len < motif_len:
            return np.array([], dtype=np.float64)

        # Map bases to indices for efficient lookup
        idx_map = {base: i for i, base in enumerate(self.alphabet)}
        
        # Prepare an array to hold the scores
        scores = np.zeros(seq_len - motif_len + 1, dtype=np.float64)

        for i in range(seq_len - motif_len + 1):
            subsequence = sequence_str[i : i + motif_len]
            score = 0.0
            for j in range(motif_len):
                base = subsequence[j]
                if base in idx_map:
                    # Lookup score from the PSSM array: row is base index, col is position
                    score += self._array[idx_map[base], j]
                # If base is not in alphabet, it's generally treated as 0 or skipped,
                # but for simplicity, we assume sequence is fully in alphabet.

            scores[i] = score

        return scores

    def search(self, sequence: Any, threshold: float = 0.0) -> List[Tuple[int, float]]:
        """Search for motif instances in a sequence that meet the score threshold."""
        
        forward_scores = self.calculate(sequence)
        results: List[Tuple[int, float]] = []

        # 1. Forward strand (positive positions)
        for pos, score in enumerate(forward_scores):
            if score >= threshold:
                results.append((pos, score))

        # 2. Reverse strand (negative positions)
        rc_pssm = self.reverse_complement()
        rc_scores = rc_pssm.calculate(sequence)
        
        # Negative positions follow Python slicing convention:
        # Position -N refers to the instance at seq[len(seq)-N : len(seq)-N + len(m)]
        # The first possible match starts at index -(seq_len - motif_len)
        for i, score in enumerate(rc_scores):
            if score >= threshold:
                # The position index for reverse strand, relative to the end of the sequence
                # i.e., position 0 of rc_scores corresponds to position -(motif_len + 0) from the end
                pos_from_end = len(sequence) - (i + self.length)
                results.append((-pos_from_end, score))
                
        return results
    
    def mean(self, background: Optional[Dict[str, float]] = None) -> float:
        """Calculate the mean score (Relative Entropy / Kullback-Leibler distance)."""
        if background is None:
            # Default uniform background
            bg_prob = 1.0 / len(self.alphabet)
            background = {base: bg_prob for base in self.alphabet}

        # Calculate KL Divergence (Relative Entropy)
        mean_score = 0.0
        # PSSM scores are log2(P_motif / P_bg)
        # E = sum_j ( sum_i ( P_motif(i,j) * PSSM(i,j) ) )
        # where P_motif is the underlying probability matrix (PWM)
        
        # Since the PSSM was calculated with an underlying PWM, we need to re-derive it
        # P_motif = background[i] * 2 ** PSSM(i,j)

        for j in range(self.length): # For each column
            column_mean = 0.0
            for i, base in enumerate(self.alphabet): # For each base
                pssm_score = self.pssm[base][j]
                bg_prob = background[base]
                
                if pssm_score == float('-inf'):
                    # If PSSM is -inf, P_motif is 0. This term contributes 0 to the sum.
                    p_motif = 0.0
                else:
                    p_motif = bg_prob * (2 ** pssm_score)
                
                # We want P_motif * PSSM_score
                if p_motif > 0 and pssm_score != float('-inf'):
                    column_mean += p_motif * pssm_score
                
            mean_score += column_mean
            
        return mean_score

    def std(self, background: Optional[Dict[str, float]] = None) -> float:
        """Calculate the standard deviation of scores."""
        if background is None:
            # Default uniform background
            bg_prob = 1.0 / len(self.alphabet)
            background = {base: bg_prob for base in self.alphabet}

        mean = self.mean(background)
        variance = 0.0

        for j in range(self.length): # For each column
            column_mean_of_squares = 0.0
            for i, base in enumerate(self.alphabet): # For each base
                pssm_score = self.pssm[base][j]
                bg_prob = background[base]
                
                if pssm_score == float('-inf'):
                    p_motif = 0.0
                else:
                    p_motif = bg_prob * (2 ** pssm_score)

                # E[X^2] = sum(p_i * x_i^2)
                if p_motif > 0 and pssm_score != float('-inf'):
                    column_mean_of_squares += p_motif * (pssm_score ** 2)

            # Var(X) = E[X^2] - E[X]^2
            # Var(sum of random variables) = sum(Var(individual variable)) if independent
            # Assuming independence between columns:
            variance += column_mean_of_squares - (mean / self.length) ** 2

        return math.sqrt(variance)

    def distribution(self, background: Optional[Dict[str, float]] = None, precision: int = 10000) -> 'Any': # Any is for ScoreDistribution
        """Calculate the score distribution for the PSSM."""
        from .thresholds import ScoreDistribution # Deferred import
        
        if background is None:
            bg_prob = 1.0 / len(self.alphabet)
            background = {base: bg_prob for base in self.alphabet}
            
        return ScoreDistribution(self, background, precision=precision)

    def reverse_complement(self) -> 'PositionSpecificScoringMatrix':
        """Return the reverse complement of the PSSM."""
        # Note: PWM provides the base implementation. We override the return type hint.
        return super().reverse_complement().log_odds(
             background={b: 0.25 for b in self.alphabet} # Log_odds implicitly uses uniform unless background is set
        )
    """