import unittest
import numpy as np
import math
import sys 
from typing import List, Dict

# Standard BioPython imports for the golden reference implementation
from Bio import motifs
from Bio.Seq import Seq

# Alias imported BioPython functions/classes to match the names used in the tests
CountsMatrix = None
PositionSpecificScoringMatrix = None
PositionWeightMatrix = None
ScoreDistribution = None

# --- DATA PATHS (Relative to this file) ---
minimal_dna_path = "../../data/minimal_test.meme"
minimal_rna_path = "../../data/minimal_test_rna.meme"

try:
    # Create a minimal motif instance to discover class types
    sequences_temp = [Seq("GATTACA"), Seq("TATTTTA")]
    motif_temp = motifs.create(sequences_temp)
    counts_temp = motif_temp.counts

    # Dynamically assign the class types
    CountsMatrix = type(counts_temp)
    
    # Get PWM type (returned by normalize)
    pwm_temp = counts_temp.normalize(pseudocounts=0.1)
    PositionWeightMatrix = type(pwm_temp)
    
    # Get PSSM type (returned by log_odds).
    pssm_temp = pwm_temp.log_odds(background={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
    PositionSpecificScoringMatrix = type(pssm_temp)
    
    # Try importing ScoreDistribution from the usual location
    from Bio.motifs import thresholds
    ScoreDistribution = thresholds.ScoreDistribution
    
except Exception as e:
    print(f"Error during dynamic class discovery (Biopython check): {e}", file=sys.stderr)
    raise e

# Re-map the names for the main motif functions
Motif = motifs.Motif
create = motifs.create
read = motifs.read


class TestMotifCreation(unittest.TestCase):
    """Tests the creation and basic properties of the Motif class using BioPython."""

    def test_create_basic(self):
        """Test creating a motif from a simple list of sequences."""
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        seq_objects = [Seq(s) for s in sequences]
        
        # Use the BioPython create function
        motif = create(seq_objects, alphabet="ACGT")

        self.assertIsInstance(motif, Motif)
        self.assertEqual(len(motif), 5)
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(motif.consensus, "AACGC")
        
        # Check counts
        self.assertAlmostEqual(motif.counts['A'][0], 3.0)
        self.assertAlmostEqual(motif.counts['C'][2], 3.0)
        self.assertAlmostEqual(motif.counts['G'][4], 2.0)
        
    def test_slice(self):
        """Test slicing the motif."""
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        seq_objects = [Seq(s) for s in sequences]
        motif = create(seq_objects, alphabet="ACGT")
        
        sub_motif = motif[1:4]
        self.assertIsInstance(sub_motif, Motif)
        self.assertEqual(len(sub_motif), 3)
        self.assertEqual(sub_motif.consensus, "ACG")
        
        # Check counts of the slice: A T A (index 1), C A G (index 2), G C C (index 3)
        # Position 1 (index 0 in sub_motif): A=3, C=0, G=0, T=1
        self.assertAlmostEqual(sub_motif.counts['A'][0], 3.0)
        self.assertAlmostEqual(sub_motif.counts['T'][0], 1.0)

    def test_reverse_complement(self):
        """Test calculating the reverse complement."""
        sequences = ["ATGC", "ATTC", "GTGC", "CTTC"]
        seq_objects = [Seq(s) for s in sequences]
        motif = create(seq_objects, alphabet="ACGT")
        
        # Consensus: ATGC
        self.assertEqual(motif.consensus, "ATGC")
        
        rc_motif = motif.reverse_complement()
        
        # RC Consensus: GCAT (Reverse complement of ATGC)
        self.assertEqual(rc_motif.consensus, "GCAT")
        self.assertEqual(len(rc_motif), 4)

        # Check counts at RC position 0 (original position 3, base G/C)
        # Original Counts at pos 3: C=4, G=0, A=0, T=0
        # RC Counts at pos 0 (C -> G, G -> C): G=4, C=0, A=0, T=0
        self.assertAlmostEqual(rc_motif.counts['G'][0], 4.0)
        self.assertAlmostEqual(rc_motif.counts['C'][0], 0.0)


class TestMatrixOperations(unittest.TestCase):
    """Tests the matrix classes (CountsMatrix, PWM, PSSM) using BioPython."""

    def setUp(self):
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        seq_objects = [Seq(s) for s in sequences]
        self.motif = create(seq_objects, alphabet="ACGT")
        self.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

    def test_counts_matrix(self):
        """Check CountsMatrix initialization and basic access."""
        counts = self.motif.counts
        # We rely on BioPython's CountsMatrix here, which is aliased above.
        self.assertIsInstance(counts, CountsMatrix) 
        self.assertEqual(counts.length, 5)
        self.assertEqual(counts.alphabet, "ACGT")
        
        self.assertAlmostEqual(counts['A'][0], 3.0)
        self.assertAlmostEqual(counts['T'][0], 1.0)
        self.assertAlmostEqual(counts['G'][4], 2.0)

    def test_pwm_creation(self):
        """Test creation of Position Weight Matrix (PWM)."""
        # Normalize counts (frequencies)
        pwm = self.motif.counts.normalize(pseudocounts=0.0)
        self.assertIsInstance(pwm, PositionWeightMatrix)
        
        # Position 0: 3 A, 1 T. Total 4.
        self.assertAlmostEqual(pwm['A'][0], 0.75)
        self.assertAlmostEqual(pwm['T'][0], 0.25)

    def test_pssm_creation(self):
        """Test creation of Position Specific Scoring Matrix (PSSM)."""
        pwm = self.motif.counts.normalize(pseudocounts=1.0) # Add pseudocounts to avoid log(0)
        pssm = pwm.log_odds(background=self.background)
        self.assertIsInstance(pssm, PositionSpecificScoringMatrix)
        
        # Position 0: 
        self.assertAlmostEqual(pssm['A'][0], 1.0)
        self.assertAlmostEqual(pssm['T'][0], 0.0)
        
    def test_pssm_score_hit(self):
        """Test scoring a sequence hit."""
        pwm = self.motif.counts.normalize(pseudocounts=1.0)
        pssm = pwm.log_odds(background=self.background)
        
        # Sequence "AACGC" (consensus)
        scores = pssm.get_scores(Seq("CAACGCGT"))
        
        self.assertAlmostEqual(scores[1], 3.17, places=2)
        
    def test_pssm_score_non_hit(self):
        """Test scoring a sequence that is not the consensus."""
        pwm = self.motif.counts.normalize(pseudocounts=1.0)
        pssm = pwm.log_odds(background=self.background)
        
        # Sequence "GATTT"
        scores = pssm.get_scores(Seq("CGATTTT"))
        self.assertAlmostEqual(scores[1], -1.0, places=2)


class TestMinimalFormatReading(unittest.TestCase):
    """Tests the minimal (MEME) format parser using BioPython."""

    def test_minimal_dna(self):
        """Test reading a DNA motif file."""
        # BioPython's read function takes an open file handle, so we wrap it.
        with open(minimal_dna_path) as handle:
            motif = read(handle, fmt="minimal")[0]
        
        self.assertIsInstance(motif, Motif)
        self.assertEqual(motif.name, "ATCCGATT") 
        self.assertEqual(len(motif), 8)
        self.assertEqual(motif.alphabet, "ACGT")
        
    def test_minimal_rna(self):
        """Test reading an RNA motif file."""
        # BioPython's read function takes an open file handle, so we wrap it.
        with open(minimal_rna_path) as handle:
            motif = read(handle, fmt="sites")[0]
        
        self.assertIsInstance(motif, Motif)
        self.assertEqual(motif.name, "KRP") 
        self.assertEqual(len(motif), 19)
        
        # Check a specific count value (Pos 0, U=10)
        self.assertAlmostEqual(motif.counts['U'][0], 10.0)


class TestScoreDistribution(unittest.TestCase):
    """Tests the ScoreDistribution class for threshold calculation using BioPython."""
    
    def setUp(self):
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        self.motif = create([Seq(s) for s in sequences], alphabet="ACGT")
        self.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        
        pwm = self.motif.counts.normalize(pseudocounts=1.0)
        self.pssm = pwm.log_odds(background=self.background)
        
        # BioPython's ScoreDistribution is used here
        self.dist = ScoreDistribution(pssm=self.pssm, background=self.background, precision=1000)

    def test_distribution_init(self):
        """Check if distribution properties are initialized."""
        self.assertTrue(self.dist.pssm.length == 5)
        # BioPython's structure uses a list of density floats
        self.assertTrue(len(self.dist.mo_density) > 0) 
        
    def test_threshold_fpr(self):
        """Test calculation of threshold based on False Positive Rate (FPR)."""
        threshold_min = self.dist.threshold_fpr(1.0)
        self.assertAlmostEqual(threshold_min, self.dist.min_score, places=2)
        
        threshold_low_fpr = self.dist.threshold_fpr(0.01)
        self.assertTrue(threshold_low_fpr > threshold_min)

    def test_threshold_balanced(self):
        """Test calculation of threshold based on balanced error rates."""
        threshold_balanced = self.dist.threshold_balanced(balance_factor=1.0)
        
        self.assertTrue(threshold_balanced > self.dist.min_score)
        self.assertTrue(threshold_balanced < self.dist.pssm.max)


if __name__ == '__main__':
    unittest.main()
