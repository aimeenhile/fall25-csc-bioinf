import unittest
import os
import sys
import numpy as np
import math
from typing import List, Dict

# DATA 
data_dir = "data" 
minimal_dna_path = data_dir + "/minimal_test.meme"
minimal_rna_path = data_dir + "/minimal_test_rna.meme"

# --- Dynamic Imports for Codon/Python Environment ---
try:
    __codon__
    CODON = True
except NameError:
    CODON = False

if CODON:
    # Codon environment
    from bio_codon import create, read, Motif
    # These are expected to be available directly within the Codon implementation namespace
    from bio_codon.matrix import CountsMatrix, PositionSpecificScoringMatrix, PositionWeightMatrix
    from bio_codon.thresholds import ScoreDistribution
    try:
        from .code.bio_codon.seq import Seq
    except ImportError:
        pass # Will rely on create to handle strings if Seq is not directly imported
else:
    # Python environment (BioPython)
    from Bio import motifs
    from Bio.Seq import Seq # Import Seq needed for object creation in setUp
    
    # FIX: The previous 'from Bio.motifs.matrix import ...' was failing.
    # We now access the required classes via the imported 'motifs' object's submodules 
    # and rename them locally to match the test case usage (e.g., CountsMatrix).
    CountsMatrix = motifs.matrix.CountsMatrix
    PositionSpecificScoringMatrix = motifs.matrix.PositionSpecificScoringMatrix
    PositionWeightMatrix = motifs.matrix.PositionWeightMatrix
    ScoreDistribution = motifs.thresholds.ScoreDistribution
    
    # Re-map the names for the main motif functions
    Motif = motifs.Motif
    create = motifs.create
    read = motifs.read


# --- Test Cases ---

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
        # A set of RNA sequences (using U instead of T)
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
        try:
            # Try creating with Seq if available (Python/Biopython)
            # Seq class is imported in Python env or expected in Codon env
            self.motif_dna = create([Seq(s) for s in self.sequences_dna], alphabet="ACGT")
            self.motif_rna = create([Seq(s) for s in self.sequences_rna], alphabet="ACGU")
        except NameError:
            # Fallback for environments where Seq is not directly imported (Codon)
            self.motif_dna = create(self.sequences_dna, alphabet="ACGT")
            self.motif_rna = create(self.sequences_rna, alphabet="ACGU")
        
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
        # This assert now works in the Python environment due to the fixed imports
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
        # Biopython breaks ties by alphabetical order (A<C<G<T). In position 3, G=3, T=3. G is chosen.
        # Expected consensus: T A C G C
        self.assertEqual(self.motif_dna.consensus, "TACGC") 
        self.assertEqual(self.motif_rna.consensus, "UACGC")

    def test_reverse_complement(self):
        """Test reverse complement property."""
        rc_motif = self.motif_dna.reverse_complement()
        self.assertIsInstance(rc_motif, Motif)
        self.assertEqual(rc_motif.length, 5)
        # Original consensus: T A C G C
        # Reverse complement: G C G T A
        self.assertEqual(rc_motif.consensus, "GCGTA")

    def test_motif_slicing(self):
        """Test slicing of the Motif object."""
        # Slice from index 1 to 4 (exclusive of 4) -> indices 1, 2, 3
        sliced_motif = self.motif_dna[1:4]
        self.assertIsInstance(sliced_motif, Motif)
        self.assertEqual(sliced_motif.length, 3)
        # Consensus of columns 1, 2, 3: A C G
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
        
        try:
            self.motif = create([Seq(s) for s in sequences], alphabet="ACGT")
        except NameError:
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
        score_t0 = math.log2((6.0 / 11.0) / 0.25) # (5+1)/11
        score_a1 = math.log2((5.0 / 11.0) / 0.25) # (4+1)/11
        score_c2 = math.log2((5.0 / 11.0) / 0.25) # (4+1)/11
        score_a3 = math.log2((1.0 / 11.0) / 0.25) # (0+1)/11
        score_a4 = math.log2((1.0 / 11.0) / 0.25) # (0+1)/11
        
        expected_score = score_t0 + score_a1 + score_c2 + score_a3 + score_a4
        
        result = pssm.calculate("TACAA")
        
        # Handle cases where calculate returns a list or a single float
        final_score = result[0] if isinstance(result, list) else result
        self.assertAlmostEqual(final_score, expected_score, places=5)

class TestMinimalParsing(unittest.TestCase):
    """Tests parsing of the minimal MEME file format."""

    def test_dna_parsing(self):
        """Test parsing of a DNA minimal MEME file."""
        with open(minimal_dna_path) as f:
            motifs_list = read(f, "minimal") 

        # We assert for 3 motifs based on the file content 
        self.assertEqual(len(motifs_list), 3) 
        
        # Check first motif (KRP)
        motif = motifs_list[0]
        self.assertEqual(motif.name, "KRP")
        self.assertEqual(motif.length, 19)
        self.assertEqual(motif.alphabet, "ACGT")

        # Check one value from the CountsMatrix
        # nsites=17, PFM[A, 4] = 0.823529. Count = round(0.823529 * 17) = 14
        self.assertAlmostEqual(motif.counts["A"][4], 14.0, delta=0.01)
        
        # Check consensus
        # Pos 0: T (0.82), Pos 1: G (0.64), Pos 2: T (0.94), Pos 3: G (0.76), Pos 4: A (0.82)
        self.assertEqual(motif.consensus[:5], "TGTGA")

    def test_rna_parsing(self):
        """Test parsing of an RNA minimal MEME file."""
        with open(minimal_rna_path) as f:
            motifs_list = read(f, "minimal") 

        # We assert for 3 motifs based on the file content 
        self.assertEqual(len(motifs_list), 3)
        
        # Check first motif (KRP)
        motif = motifs_list[0]
        self.assertEqual(motif.name, "KRP") 
        self.assertEqual(motif.length, 19)
        self.assertEqual(motif.alphabet, "ACGT") # MEME files use T even for RNA, Biopython handles the conversion
        
        # Check one value from the CountsMatrix (T is counted as U internally by Biopython for RNA)
        # nsites=17, PFM[T, 0] = 0.823529. Count = round(0.823529 * 17) = 14
        self.assertAlmostEqual(motif.counts["T"][0], 14.0, delta=0.01)

        # Check consensus - Biopython's read function usually leaves the alphabet as ACGT 
        # unless told otherwise, but it uses the T/U conversion internally for consensus
        # Pos 0: T (0.82), Pos 1: G (0.64), Pos 2: T (0.94), Pos 3: G (0.76), Pos 4: A (0.82)
        self.assertEqual(motif.consensus[:5], "TGTGA")


class TestScoreDistribution(unittest.TestCase):
    """Tests the ScoreDistribution class for threshold calculation."""
    
    def setUp(self):
        sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
        
        try:
            self.motif = create([Seq(s) for s in sequences], alphabet="ACGT")
        except NameError:
            self.motif = create(sequences, alphabet="ACGT")
            
        self.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        
        # PSSM is required for ScoreDistribution
        pwm = self.motif.counts.normalize(pseudocounts=1.0)
        self.pssm = pwm.log_odds(background=self.background)
        
        # Initialize ScoreDistribution
        self.dist = ScoreDistribution(pssm=self.pssm, background=self.background, precision=10000)

    def test_distribution_init(self):
        """Check if distribution properties are initialized."""
        # Check some basic properties from the initialization
        self.assertTrue(self.dist.interval > 0.0)
        self.assertTrue(self.dist.n_points > 0)
        self.assertTrue(self.dist.ic is not None) 
        # Check that the densities are calculated
        self.assertTrue(len(self.dist.mo_density) > 0)
        self.assertTrue(len(self.dist.bg_density) > 0)
        
    def test_threshold_fpr(self):
        """Test calculation of threshold based on False Positive Rate (FPR)."""
        # A low FPR should result in a high score threshold
        threshold_low_fpr = self.dist.threshold_fpr(0.01)
        # A high FPR should result in a low score threshold
        threshold_high_fpr = self.dist.threshold_fpr(0.5)

        self.assertTrue(threshold_low_fpr > threshold_high_fpr)
        # Ensure the threshold is within the possible range (min_score to max_score of the PSSM).
        self.assertTrue(self.pssm.min < threshold_low_fpr < self.pssm.max)
        self.assertTrue(self.pssm.min < threshold_high_fpr < self.pssm.max)


if __name__ == '__main__':
    # Use unittest.main() to run all tests
    unittest.main(argv=sys.argv[:1], exit=False)