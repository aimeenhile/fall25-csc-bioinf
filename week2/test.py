import unittest
import os
import sys
import numpy as np
import math
from typing import List, Dict, TextIO
from Bio.Seq import Seq

DATA = "data"
data_dir = os.path.join(os.path.dirname(__file__), "data")
minimal_dna_path = os.path.join(data_dir, "minimal_test.meme")
minimal_rna_path = os.path.join(data_dir, "minimal_test_rna.meme")

"""
try:
    __codon__
    CODON = True
except NameError:
    CODON = False

if CODON:
    # Codon environment
    from bio_codon import create, read, Motif
    from bio_codon.matrix import PositionSpecificScoringMatrix, PositionWeightMatrix, CountsMatrix
    from bio_codon.thresholds import ScoreDistribution
    try:
        from .code.bio_codon.seq import Seq
    except ImportError:
        pass # Will rely on create to handle strings if Seq is not directly imported
else:
    # Python environment
    from Bio import motifs
    from Bio.motifs import matrix, thresholds 

    # Re-map the names to match the variables used in the test cases for consistency
    Motif = motifs.Motif
    create = motifs.create
    read = motifs.read
"""

    
if hasattr(str, 'memcpy'):
    # --- Codon Environment Imports ---
    # When running 'codon run test.py' from week2/, imports must reference the 
    # top-level directory ('bio_codon').
    from bio_codon.motifs import create, read, Motif
    from bio_codon.motifs.matrix import PositionSpecificScoringMatrix, PositionWeightMatrix, CountsMatrix
    from bio_codon.motifs.thresholds import ScoreDistribution
    
    # The global 'from Bio.Seq import Seq' should suffice if the Codon bridge is working.
    
else:
    # --- Python Environment Imports ---
    # We use the standard, installed Biopython library for reference tests.
    from Bio import motifs
    # Import matrix/threshold classes directly from their explicit submodules
    import Bio.motifs.matrix
    import Bio.motifs.thresholds
    # Re-map the names to match the variables used in the test cases for consistency
    Motif = motifs.Motif
    create = motifs.create
    read = motifs.read

    CountsMatrix = Bio.motifs.matrix.CountsMatrix
    PositionWeightMatrix = Bio.motifs.matrix.PositionWeightMatrix
    PositionSpecificScoringMatrix = Bio.motifs.matrix.PositionSpecificScoringMatrix
    ScoreDistribution = Bio.motifs.thresholds.ScoreDistribution

    
# --- Test Cases ---

import math
import tempfile
import unittest

try:
    import numpy as np
except ImportError:
    from Bio import MissingExternalDependencyError

    raise MissingExternalDependencyError(
        "Install numpy if you want to use Bio.motifs."
    ) from None

from Bio import motifs
from Bio.Seq import Seq


class TestMotifBasic(unittest.TestCase):
    """Tests basic Motif creation and properties."""
    
    def setUp(self):
        # A set of DNA sequences to create a motif
        self.sequences_dna = [
            "TACAA",
            "TACGC",
            "TACAC",
            "TACCC",
            "AACCC",
            "AATGC",
            "AATGC",
        ]
        # A set of RNA sequences
        self.sequences_rna = [
            "UACAA",
            "UACGC",
            "UACAC",
            "UACCC",
            "AACCC",
            "AAUGC",
            "AAUGC",
        ]
        
        # Create Motif objects using the custom `create` function
        self.motif_dna = create([Seq(s) for s in self.sequences_dna], alphabet="ACGT")
        self.motif_rna = create([Seq(s) for s in self.sequences_rna], alphabet="ACGU")
        
        # Expected Counts Matrix for DNA motif (Position: 0 1 2 3 4)
        # A: 2 4 0 0 0
        # C: 0 0 4 1 4
        # G: 0 0 0 3 0
        # T: 5 3 3 3 3 (T/U)
        
    def test_motif_length(self):
        """Check if motif length is calculated correctly."""
        self.assertEqual(self.motif_dna.length, 5)
        self.assertEqual(self.motif_rna.length, 5)

    def test_motif_count_matrix(self):
        """Check if the CountsMatrix attribute is correct."""
        counts = self.motif_dna.counts
        self.assertIsInstance(counts, CountsMatrix)
        self.assertEqual(counts["A", 0], 2)
        self.assertEqual(counts["T", 0], 5)
        self.assertEqual(counts["C", 2], 4)
        self.assertEqual(counts["G", 3], 3)
        self.assertEqual(counts["A", 4], 0)
        
        counts_rna = self.motif_rna.counts
        # For RNA, U replaces T
        self.assertEqual(counts_rna["U", 0], 5)

    def test_motif_consensus(self):
        """Check if consensus sequence is generated correctly."""
        # T/A (T dominant), A/T/C (A dominant), C/T (C dominant), G/C/T (G/T dominant, but G is 3, T is 3)
        # A: 2 4 0 0 0
        # C: 0 0 4 1 4
        # G: 0 0 0 3 0
        # T: 5 3 3 3 3
        # T A C G T -> T A C T C (If ties are broken by alphabet order, G comes before T)
        # Biopython breaks ties by alphabetical order (A<C<G<T). In position 3, G=3, T=3. G is chosen.
        # In position 4, C=4, T=3. C is chosen.
        # The expected Biopython consensus (which should be replicated):
        # Pos 0: T (5) > A (2)
        # Pos 1: A (4) > T (3)
        # Pos 2: C (4) > T (3)
        # Pos 3: G (3) = T (3). Tie-break: G is before T. G chosen.
        # Pos 4: C (4) > T (3)
        self.assertEqual(self.motif_dna.consensus, "T A C G C".replace(" ", "")) 
        self.assertEqual(self.motif_rna.consensus, "U A C G C".replace(" ", ""))

    def test_reverse_complement(self):
        """Test reverse complement property."""
        rc_motif = self.motif_dna.reverse_complement()
        self.assertIsInstance(rc_motif, Motif)
        self.assertEqual(rc_motif.length, 5)
        # Original consensus: T A C G C
        # Reverse complement: G C G T A
        self.assertEqual(rc_motif.consensus, "G C G T A".replace(" ", ""))

    def test_motif_slicing(self):
        """Test slicing of the Motif object."""
        # Slice from index 1 to 4 (exclusive of 4) -> indices 1, 2, 3
        sliced_motif = self.motif_dna[1:4]
        self.assertIsInstance(sliced_motif, Motif)
        self.assertEqual(sliced_motif.length, 3)
        # Original columns 1, 2, 3:
        # T A C G C
        # A: 4 0 0
        # C: 0 4 1
        # G: 0 0 3
        # T: 3 3 3
        # Consensus: A C G
        self.assertEqual(sliced_motif.consensus, "ACG")

class TestMotifMatrixConversions(unittest.TestCase):
    """Tests matrix conversions (PWM, PSSM) and scoring."""
    
    def setUp(self):
        # Use the same sequences as TestMotifBasic
        sequences = [
            "TACAA",
            "TACGC",
            "TACAC",
            "TACCC",
            "AACCC",
            "AATGC",
            "AATGC",
        ]
        self.motif = create(sequences, alphabet="ACGT")
        # Standard Biopython default background is equiprobable (0.25 for each)
        self.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

    def test_pwm_conversion(self):
        """Test conversion to Position Weight Matrix (PWM)."""
        # normalize() returns a PositionWeightMatrix, which is the Frequency Matrix in Biopython terms
        pwm = self.motif.counts.normalize(pseudocounts=1)
        self.assertIsInstance(pwm, PositionWeightMatrix)
        self.assertEqual(pwm.length, 5)
        
        # Total counts at pos 0 is 7. With pseudocount=1, total is 7+4=11.
        # A count is 2. Normalized freq = (2+1)/11 = 3/11
        self.assertAlmostEqual(pwm["A", 0], 3/11, places=5)
        # T count is 5. Normalized freq = (5+1)/11 = 6/11
        self.assertAlmostEqual(pwm["T", 0], 6/11, places=5)
        
    def test_pssm_conversion(self):
        """Test conversion to Position Specific Scoring Matrix (PSSM)."""
        # PSSM is log_odds(PWM)
        pwm = self.motif.counts.normalize(pseudocounts=1)
        pssm = pwm.log_odds(background=self.background)
        self.assertIsInstance(pssm, PositionSpecificScoringMatrix)
        
        # Check PSSM calculation: PSSM[A, 0] = log2( (2+1)/(7+4) / 0.25 ) = log2( (3/11) / 0.25 )
        expected_score_a0 = math.log2((3.0 / 11.0) / 0.25)
        self.assertAlmostEqual(pssm["A", 0], expected_score_a0, places=5)
        
        # Check PSSM calculation: PSSM[T, 0] = log2( (5+1)/(7+4) / 0.25 ) = log2( (6/11) / 0.25 )
        expected_score_t0 = math.log2((6.0 / 11.0) / 0.25)
        self.assertAlmostEqual(pssm["T", 0], expected_score_t0, places=5)

    def test_pssm_calculate_score(self):
        """Test PSSM scoring of a sequence."""
        pwm = self.motif.counts.normalize(pseudocounts=1)
        pssm = pwm.log_odds(background=self.background)
        
        # Score the first sequence "TACAA" (length 5)
        # Expected score: PSSM[T,0] + PSSM[A,1] + PSSM[C,2] + PSSM[A,3] + PSSM[A,4]
        score_t0 = math.log2((6.0 / 11.0) / 0.25)
        score_a1 = math.log2((5.0 / 11.0) / 0.25)
        score_c2 = math.log2((5.0 / 11.0) / 0.25)
        score_a3 = math.log2((1.0 / 11.0) / 0.25) # A count 0 -> (0+1)/11
        score_a4 = math.log2((1.0 / 11.0) / 0.25) # A count 0 -> (0+1)/11
        
        expected_score = score_t0 + score_a1 + score_c2 + score_a3 + score_a4
        
        result = pssm.calculate("TACAA")
        # In Biopython's PSSM, 'calculate' returns a list of scores for all possible start positions.
        # Since the motif length is 5 and sequence length is 5, there is only one start position (index 0).
        if isinstance(result, list) and len(result) == 1:
             self.assertAlmostEqual(result[0], expected_score, places=5)
        elif isinstance(result, float):
             # Handle a single float return if the implementation simplifies the return for length match
             self.assertAlmostEqual(result, expected_score, places=5)
        else:
             self.assertAlmostEqual(result[0], expected_score, places=5)

class TestMinimalParsing(unittest.TestCase):
    """Tests parsing of the minimal MEME file format."""

    def test_dna_parsing(self):
        """Test parsing of a DNA minimal MEME file."""
        if not os.path.exists(minimal_dna_path):
            print(f"Warning: DNA test file not found at {minimal_dna_path}")
            return
            
        with open(minimal_dna_path) as handle:
            # Use the local `read` function
            motifs_list = read(handle, "minimal") 

        self.assertEqual(len(motifs_list), 2)
        
        motif = motifs_list[0]
        self.assertIsInstance(motif, Motif)
        self.assertEqual(motif.name, "KRP")
        self.assertEqual(motif.length, 19)
        self.assertAlmostEqual(motif.evalue, 4.1e-09, delta=1e-11)
        
        # Check consensus of the first motif (from minimal_test.meme)
        # Pos 0: T (0.82)
        # Pos 1: G (0.64)
        # Pos 2: T (0.94)
        # Pos 3: G (0.76)
        # Pos 4: A (0.82)
        self.assertEqual(motif.consensus[0:5], "TGTGA")

        # Check second motif
        motif_2 = motifs_list[1]
        self.assertEqual(motif_2.name, "IFXA")
        self.assertEqual(motif_2.length, 25)

    def test_rna_parsing(self):
        """Test parsing of an RNA minimal MEME file."""
        if not os.path.exists(minimal_rna_path):
            print(f"Warning: RNA test file not found at {minimal_rna_path}")
            return
            
        with open(minimal_rna_path) as handle:
            # Use the local `read` function
            motifs_list = read(handle, "minimal") 

        self.assertEqual(len(motifs_list), 2)
        
        motif = motifs_list[0]
        self.assertIsInstance(motif, Motif)
        self.assertEqual(motif.name, "KRP_fake_RNA")
        self.assertEqual(motif.length, 19)
        # Check consensus - should use U instead of T
        # Pos 0: U (0.82)
        # Pos 1: G (0.64)
        # Pos 2: U (0.94)
        # Pos 3: G (0.76)
        # Pos 4: A (0.82)
        self.assertEqual(motif.consensus[0:5], "UGUGA")


class TestScoreDistribution(unittest.TestCase):
    """Tests the ScoreDistribution class for threshold calculation."""
    
    def setUp(self):
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        self.motif = create(sequences, alphabet="ACGT")
        self.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        
        # PSSM is required for ScoreDistribution
        pwm = self.motif.counts.normalize(pseudocounts=1.0)
        self.pssm = pwm.log_odds(background=self.background)
        
        # Initialize ScoreDistribution
        self.dist = ScoreDistribution(pssm=self.pssm, background=self.background, precision=1000)

    def test_distribution_init(self):
        """Check if distribution properties are initialized."""
        # Check some basic properties from the initialization
        self.assertTrue(self.dist.interval > 0.0)
        self.assertTrue(self.dist.n_points > 0)
        self.assertTrue(self.dist.ic is not None)
        self.assertEqual(len(self.dist.mo_density), self.dist.n_points)
        self.assertEqual(len(self.dist.bg_density), self.dist.n_points)
        
    def test_threshold_fpr(self):
        """Test calculation of threshold based on False Positive Rate (FPR)."""
        # A low FPR should result in a high score threshold
        threshold_low_fpr = self.dist.threshold_fpr(0.01)
        # A high FPR should result in a low score threshold
        threshold_high_fpr = self.dist.threshold_fpr(0.5)

        self.assertTrue(threshold_low_fpr > threshold_high_fpr)
        # Note: We can't check the exact value without replicating the dynamic programming,
        # but we can ensure the threshold is within the possible range (min_score to max_score of the PSSM).
        self.assertTrue(self.pssm.min < threshold_low_fpr < self.pssm.max)
        self.assertTrue(self.pssm.min < threshold_high_fpr < self.pssm.max)


if __name__ == '__main__':
    # Use unittest.main() to run all tests
    unittest.main(argv=sys.argv[:1], exit=False)