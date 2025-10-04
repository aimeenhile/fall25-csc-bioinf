import unittest
import numpy as np
import math
from typing import List, Dict

# --- DATA PATHS (relative to week2/bio_codon/motifs/) ---
# Go up two levels (../../) to reach the 'week2/' root, then down into 'data/'
minimal_dna_path = "../../data/minimal_test.meme"
minimal_rna_path = "../../data/minimal_test_rna.meme"

# --- Codon Environment Imports ---
# Assumes this test file is run within the package structure
from . import create, read, Motif 
from ..matrix import CountsMatrix, PositionSpecificScoringMatrix, PositionWeightMatrix
from ..thresholds import ScoreDistribution

# Import the Codon-wrapped version of BioPython's Seq
from python import Bio as cBio
Seq = cBio.Seq.Seq


class TestMotifBasics(unittest.TestCase):
    """Tests the core Motif class methods (consensus, RC, slicing)."""
    
    def setUp(self):
        # Create a simple motif for testing
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        self.seqs = [Seq(s) for s in sequences]
        self.motif = create(self.seqs, alphabet="ACGT")

    def test_motif_creation(self):
        """Check if motif object is correctly initialized."""
        self.assertEqual(len(self.motif), 5)
        self.assertEqual(self.motif.length, 5)
        self.assertEqual(self.motif.alphabet, "ACGT")
        self.assertEqual(len(self.motif.instances), 4)

    def test_consensus(self):
        """Test consensus sequence calculation."""
        # Position counts:
        # A: 3, 1, 0, 0, 0
        # C: 1, 0, 3, 0, 1
        # G: 0, 0, 1, 1, 2
        # T: 0, 3, 0, 3, 1
        # Consensus: T, A, C/A, T, G/C -> should resolve ties alphabetically/by input order 
        # For our simple code, 'T', 'A', 'C', 'T', 'G' (C wins tie at pos 2, G wins tie at pos 4)
        self.assertEqual(self.motif.consensus, "TACCGC") # Recalculating based on counts:
        # Pos 0: T(1) A(3) C(0) G(0) -> A
        # Pos 1: A(1) C(0) G(0) T(3) -> T
        # Pos 2: A(0) C(3) G(1) T(0) -> C
        # Pos 3: A(0) C(0) G(1) T(3) -> T
        # Pos 4: A(0) C(1) G(2) T(1) -> G
        self.assertEqual(self.motif.consensus, "ATCTG")
        
    def test_reverse_complement(self):
        """Test reverse complement calculation."""
        rc_motif = self.motif.reverse_complement()
        self.assertEqual(len(rc_motif), 5)
        self.assertEqual(rc_motif.consensus, "CAGAT") # RC of "ATCTG" is "CAGAT"


class TestMatrix(unittest.TestCase):
    """Tests the matrix classes (CountsMatrix, PWM, PSSM)."""
    
    def setUp(self):
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        self.motif = create([Seq(s) for s in sequences], alphabet="ACGT")
        self.counts = self.motif.counts
        self.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

    def test_counts_matrix(self):
        """Check basic counts and indexing."""
        self.assertEqual(self.counts['A', 0], 3.0)
        self.assertEqual(self.counts['T', 1], 3.0)
        self.assertEqual(self.counts['C', 2], 3.0)
        self.assertEqual(self.counts['G', 4], 2.0)
        self.assertEqual(len(self.counts), 5)

    def test_normalize_to_pwm(self):
        """Test normalization with pseudocounts."""
        # Total counts per column: 4, 4, 4, 4, 4
        # Add 1.0 pseudocount to each: Total becomes 8 per column.
        pwm = self.counts.normalize(pseudocounts=1.0)
        
        # Original 'A' count at pos 0 was 3.0. With PC 1.0: (3.0 + 1.0) / 8.0 = 0.5
        self.assertAlmostEqual(pwm['A', 0], 0.5)
        
        # Original 'T' count at pos 2 was 0.0. With PC 1.0: (0.0 + 1.0) / 8.0 = 0.125
        self.assertAlmostEqual(pwm['T', 2], 0.125)
        
        # Check that one row sums to 1.0 (it shouldn't, rows are not probabilities)
        # Check that one column sums to 1.0 (it should)
        col_sum = sum(pwm[base, 1] for base in "ACGT")
        self.assertAlmostEqual(col_sum, 1.0)

    def test_pwm_to_pssm(self):
        """Test log-odds calculation."""
        pwm = self.counts.normalize(pseudocounts=1.0)
        pssm = pwm.log_odds(background=self.background)
        
        # PWM['A', 0] = 0.5. Background['A'] = 0.25.
        # Score = log2(0.5 / 0.25) = log2(2) = 1.0
        self.assertAlmostEqual(pssm['A', 0], 1.0)

        # PWM['T', 2] = 0.125. Background['T'] = 0.25.
        # Score = log2(0.125 / 0.25) = log2(0.5) = -1.0
        self.assertAlmostEqual(pssm['T', 2], -1.0)
        
    def test_pssm_scoring(self):
        """Test scoring of sequences using PSSM."""
        pwm = self.counts.normalize(pseudocounts=1.0)
        pssm = pwm.log_odds(background=self.background)
        
        # Sequence: "TATAC" (Length 5)
        # T @ 0: P(T,0)=1.0, A @ 1: P(A,1)=-0.415, T @ 2: P(T,2)=-1.0, A @ 3: P(A,3)=-0.415, C @ 4: P(C,4)=-0.415
        # Total score for TATAC should be 1.0 - 0.415 - 1.0 - 0.415 - 0.415 = -1.245 (approx, based on 1.0 PC)
        # Let's use a known high scoring sequence: "ATCTG" (the consensus)
        # Scores:
        # A@0: 1.0
        # T@1: 1.0
        # C@2: 1.0
        # T@3: 1.0
        # G@4: 0.415 (G=3+1=4/8=0.5, log2(0.5/0.25)=1.0... Wait, G count is 2. (2+1)/8 = 0.375. log2(0.375/0.25) = 0.585)
        # G@4 score should be log2(((2.0+1.0)/8.0) / 0.25) = log2(0.375/0.25) = log2(1.5) approx 0.585
        # ATCTG Score (using approximate values): 1.0 + 1.0 + 1.0 + 1.0 + 0.585 = 4.585
        
        high_score_seq = "ATCTG" # The consensus sequence
        scores = pssm.calculate(high_score_seq)
        
        self.assertEqual(len(scores), 1)
        self.assertAlmostEqual(scores[0], 4.5849625, places=5)
        
        # Test unknown base (should result in NaN)
        scores_nan = pssm.calculate("TAXCG")
        self.assertEqual(len(scores_nan), 1)
        self.assertTrue(math.isnan(scores_nan[0]))

    def test_pssm_sliding_window(self):
        """Test PSSM scoring across a longer sequence."""
        # Motif length is 5
        long_seq = "CATCTGA" # 7 bases, two windows of length 5: CATCT, ATCTG, TCTGA
        
        pwm = self.counts.normalize(pseudocounts=1.0)
        pssm = pwm.log_odds(background=self.background)
        
        scores = pssm.calculate(long_seq)
        self.assertEqual(len(scores), 3)
        
        # CATCT (Expected: 1.0 + (-0.415) + 1.0 + 1.0 + (-1.0) = 1.585)
        self.assertAlmostEqual(scores[0], 1.5849625, places=5)
        # ATCTG (Expected: 4.585)
        self.assertAlmostEqual(scores[1], 4.5849625, places=5)
        
        
class TestMinimalParser(unittest.TestCase):
    """Tests the parsing of the MEME minimal format."""
    
    def test_parse_minimal_dna(self):
        """Test parsing a single DNA motif from file."""
        # Note: The 'read' function is expected to return a single Motif object
        motif = read(minimal_dna_path, "minimal")
        
        # Check basic properties
        self.assertIsInstance(motif, Motif)
        self.assertEqual(motif.name, "KRP") 
        self.assertEqual(motif.length, 19)
        self.assertEqual(motif.alphabet, "ACGT")
        
        # Check consensus (TGTGA... from the data file)
        self.assertEqual(motif.consensus[:5], "TGTGA")

    def test_parse_minimal_rna(self):
        """Test parsing a single RNA motif from file."""
        motif = read(minimal_rna_path, "minimal")
        
        self.assertIsInstance(motif, Motif)
        self.assertEqual(motif.name, "RNA_MOTIF") 
        self.assertEqual(motif.length, 6)
        # Alphabet is ACGU based on the file content
        self.assertEqual(motif.alphabet, "ACGU") 
        
        # Check consensus (AUCCGG)
        self.assertEqual(motif.consensus, "AUCCGG")


class TestScoreDistribution(unittest.TestCase):
    """Tests the ScoreDistribution class for threshold calculation."""
    
    def setUp(self):
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        self.motif = create([Seq(s) for s in sequences], alphabet="ACGT")
        self.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        
        pwm = self.motif.counts.normalize(pseudocounts=1.0)
        self.pssm = pwm.log_odds(background=self.background)
        
        # Use low precision for faster test execution if needed, but 10k is standard
        self.dist = ScoreDistribution(pssm=self.pssm, background=self.background, precision=10000)

    def test_distribution_init(self):
        """Check if distribution properties are initialized."""
        # Check that the total probability is close to 1.0
        total_prob = sum(self.dist.distribution.values())
        self.assertAlmostEqual(total_prob, 1.0, places=4)
        
    def test_threshold_fpr(self):
        """Test calculation of threshold based on False Positive Rate (FPR)."""
        # The threshold should be high for a low FPR
        threshold_low_fpr = self.dist.threshold_fpr(0.001)
        # The threshold should be low for a high FPR
        threshold_high_fpr = self.dist.threshold_fpr(0.99)
        
        self.assertTrue(threshold_low_fpr > threshold_high_fpr)
        # Check that the max score is less than the theoretical max score (4.585)
        self.assertTrue(threshold_low_fpr <= 4.585) 

    def test_threshold_fnr(self):
        """Test calculation of threshold based on False Negative Rate (FNR)."""
        # The threshold should be low for a low FNR (most true motifs score above it)
        threshold_low_fnr = self.dist.threshold_fnr(0.001)
        # The threshold should be high for a high FNR (most true motifs score below it)
        threshold_high_fnr = self.dist.threshold_fnr(0.99)
        
        self.assertTrue(threshold_low_fnr < threshold_high_fnr)
        
    def test_threshold_balanced(self):
        """Test calculation of a balanced threshold (FPR ~ FNR)."""
        # The balanced threshold should fall somewhere between the min and max scores
        balanced_threshold = self.dist.threshold_balanced(balance_factor=1.0)
        
        self.assertTrue(balanced_threshold > self.dist.pssm.min)
        self.assertTrue(balanced_threshold < self.dist.pssm.max)


if __name__ == '__main__':
    # When running under Codon, the test discovery mechanism might handle this, 
    # but include for direct execution support.
    unittest.main()
