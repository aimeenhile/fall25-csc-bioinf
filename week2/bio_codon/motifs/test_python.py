import unittest
import numpy as np
import math
import sys 
from typing import List, Dict

# --- DATA PATHS (relative to week2/bio_codon/motifs/) ---
# Go up two levels (../../) to reach the 'week2/' root, then down into 'data/'
minimal_dna_path = "../../data/minimal_test.meme"
minimal_rna_path = "../../data/minimal_test_rna.meme"

# --- Python Environment Imports (Standard package structure) ---
# When run as part of a package, these relative imports should work in Python 3
try:
    from . import create, read, Motif 
    from ..matrix import CountsMatrix, PositionSpecificScoringMatrix, PositionWeightMatrix
    from ..thresholds import ScoreDistribution
except ImportError:
    # Fallback for direct execution of this script outside a proper package setup
    sys.path.append('..')
    from motifs import create, read, Motif
    from motifs.matrix import CountsMatrix, PositionSpecificScoringMatrix, PositionWeightMatrix
    from motifs.thresholds import ScoreDistribution

# Use the standard Biopython Seq class
from Bio.Seq import Seq


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
        # Recalculated consensus based on counts
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
        
        self.assertAlmostEqual(pwm['A', 0], 0.5)
        self.assertAlmostEqual(pwm['T', 2], 0.125)
        
        col_sum = sum(pwm[base, 1] for base in "ACGT")
        self.assertAlmostEqual(col_sum, 1.0)

    def test_pwm_to_pssm(self):
        """Test log-odds calculation."""
        pwm = self.counts.normalize(pseudocounts=1.0)
        pssm = pwm.log_odds(background=self.background)
        
        # log2(0.5 / 0.25) = 1.0
        self.assertAlmostEqual(pssm['A', 0], 1.0)

        # log2(0.125 / 0.25) = -1.0
        self.assertAlmostEqual(pssm['T', 2], -1.0)
        
    def test_pssm_scoring(self):
        """Test scoring of sequences using PSSM."""
        pwm = self.counts.normalize(pseudocounts=1.0)
        pssm = pwm.log_odds(background=self.background)
        
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
        long_seq = "CATCTGA"
        
        pwm = self.counts.normalize(pseudocounts=1.0)
        pssm = pwm.log_odds(background=self.background)
        
        scores = pssm.calculate(long_seq)
        self.assertEqual(len(scores), 3)
        
        self.assertAlmostEqual(scores[0], 1.5849625, places=5)
        self.assertAlmostEqual(scores[1], 4.5849625, places=5)
        
        
class TestMinimalParser(unittest.TestCase):
    """Tests the parsing of the MEME minimal format."""
    
    def test_parse_minimal_dna(self):
        """Test parsing a single DNA motif from file."""
        motif = read(minimal_dna_path, "minimal")
        
        self.assertIsInstance(motif, Motif)
        self.assertEqual(motif.name, "KRP") 
        self.assertEqual(motif.length, 19)
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(motif.consensus[:5], "TGTGA")

    def test_parse_minimal_rna(self):
        """Test parsing a single RNA motif from file."""
        motif = read(minimal_rna_path, "minimal")
        
        self.assertIsInstance(motif, Motif)
        self.assertEqual(motif.name, "RNA_MOTIF") 
        self.assertEqual(motif.length, 6)
        self.assertEqual(motif.alphabet, "ACGU") 
        self.assertEqual(motif.consensus, "AUCCGG")


class TestScoreDistribution(unittest.TestCase):
    """Tests the ScoreDistribution class for threshold calculation."""
    
    def setUp(self):
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        self.motif = create([Seq(s) for s in sequences], alphabet="ACGT")
        self.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        
        pwm = self.motif.counts.normalize(pseudocounts=1.0)
        self.pssm = pwm.log_odds(background=self.background)
        
        self.dist = ScoreDistribution(pssm=self.pssm, background=self.background, precision=10000)

    def test_distribution_init(self):
        """Check if distribution properties are initialized."""
        total_prob = sum(self.dist.distribution.values())
        self.assertAlmostEqual(total_prob, 1.0, places=4)
        
    def test_threshold_fpr(self):
        """Test calculation of threshold based on False Positive Rate (FPR)."""
        threshold_low_fpr = self.dist.threshold_fpr(0.001)
        threshold_high_fpr = self.dist.threshold_fpr(0.99)
        
        self.assertTrue(threshold_low_fpr > threshold_high_fpr)
        self.assertTrue(threshold_low_fpr <= 4.585) 

    def test_threshold_fnr(self):
        """Test calculation of threshold based on False Negative Rate (FNR)."""
        threshold_low_fnr = self.dist.threshold_fnr(0.001)
        threshold_high_fnr = self.dist.threshold_fnr(0.99)
        
        self.assertTrue(threshold_low_fnr < threshold_high_fnr)
        
    def test_threshold_balanced(self):
        """Test calculation of a balanced threshold (FPR ~ FNR)."""
        balanced_threshold = self.dist.threshold_balanced(balance_factor=1.0)
        
        self.assertTrue(balanced_threshold > self.dist.pssm.min)
        self.assertTrue(balanced_threshold < self.dist.pssm.max)


if __name__ == '__main__':
    # When running directly with Python
    unittest.main()
