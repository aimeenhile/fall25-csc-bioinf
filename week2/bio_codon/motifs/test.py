import unittest
import numpy as np
import math
import sys 
from typing import List, Dict

# --- DATA PATHS (Corrected for new location: bio_codon/motifs/test.py) ---
# We must go up two levels (../../) to reach the 'week2/' root, then down into 'data/'
minimal_dna_path = "../../data/minimal_test.meme"
minimal_rna_path = "../../data/minimal_test_rna.meme"

# --- Dynamic Imports for Codon/Python Environment ---
try:
    __codon__
    CODON = True
except NameError:
    CODON = False

if CODON:
    # Codon environment: Use relative imports to find modules within the 'bio_codon' package.
    try:
        # Assuming create, read, Motif are exposed at the bio_codon package root (../)
        from .. import create, read, Motif 
        # Assuming matrix and thresholds are sibling directories to motifs
        from ..matrix import CountsMatrix, PositionSpecificScoringMatrix, PositionWeightMatrix
        from ..thresholds import ScoreDistribution
        from ..seq import Seq 
    except ImportError as e:
        print(f"Codon import error. Check that create, read, Motif, etc., are correctly exposed in the 'bio_codon' package: {e}", file=sys.stderr)
        raise e
else:
    # Python environment (BioPython): Dynamic class discovery for comparison
    from Bio import motifs
    from Bio.Seq import Seq
    
    # --- DYNAMIC CLASS DISCOVERY (Biopython) ---
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
        
        # Create Motif objects using the custom `create` function, always using Seq objects
        self.motif_dna = create([Seq(s) for s in self.sequences_dna], alphabet="ACGT")
        self.motif_rna = create([Seq(s) for s in self.sequences_rna], alphabet="ACGU")
        
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
        # Check T/U handling based on the environment
        if CODON and "U" in counts_rna.alphabet:
            self.assertEqual(counts_rna["U", 0], 5)
        elif not CODON:
            self.assertEqual(counts_rna["U", 0], 5)


    def test_motif_consensus(self):
        """Check if consensus sequence is generated correctly."""
        # Expected consensus: T A C G C
        self.assertEqual(self.motif_dna.consensus, "TACGC") 
        # Expected consensus: U A C G C
        self.assertEqual(self.motif_rna.consensus, "UACGC")

    def test_reverse_complement(self):
        """Test reverse complement property."""
        rc_motif = self.motif_dna.reverse_complement()
        self.assertIsInstance(rc_motif, Motif)
        self.assertEqual(rc_motif.length, 5)
        # Original consensus: T A C G C -> Reverse complement: G C G T A
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
        sequences = [
            "TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC",
        ]
        self.motif = create([Seq(s) for s in sequences], alphabet="ACGT")
        self.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

    def test_pwm_conversion(self):
        """Test conversion to Position Weight Matrix (PWM)."""
        pwm = self.motif.counts.normalize(pseudocounts=1)
        self.assertIsInstance(pwm, PositionWeightMatrix)
        # T count is 5. Total counts is 7. Total normalized is 7+4=11. Normalized freq = (5+1)/11 = 6/11
        self.assertAlmostEqual(pwm["T", 0], 6/11, places=5)
        
    def test_pssm_conversion(self):
        """Test conversion to Position Specific Scoring Matrix (PSSM)."""
        pwm = self.motif.counts.normalize(pseudocounts=1)
        pssm = pwm.log_odds(background=self.background)
        self.assertIsInstance(pssm, PositionSpecificScoringMatrix)
        
        # Check PSSM calculation: PSSM[T, 0] = log2( (6/11) / 0.25 )
        expected_score_t0 = math.log2((6.0 / 11.0) / 0.25)
        self.assertAlmostEqual(pssm["T", 0], expected_score_t0, places=5)

    def test_pssm_calculate_score(self):
        """Test PSSM scoring of a sequence."""
        pwm = self.motif.counts.normalize(pseudocounts=1)
        pssm = pwm.log_odds(background=self.background)
        
        # Calculate expected score for "TACAA"
        score_t0 = math.log2((6.0 / 11.0) / 0.25) 
        score_a1 = math.log2((5.0 / 11.0) / 0.25) 
        score_c2 = math.log2((5.0 / 11.0) / 0.25) 
        score_a3 = math.log2((1.0 / 11.0) / 0.25) 
        score_a4 = math.log2((1.0 / 11.0) / 0.25) 
        expected_score = score_t0 + score_a1 + score_c2 + score_a3 + score_a4
        
        result = pssm.calculate("TACAA")
        final_score = result[0] if isinstance(result, list) else result
        self.assertAlmostEqual(final_score, expected_score, places=5)

class TestMinimalParsing(unittest.TestCase):
    """Tests parsing of the minimal MEME file format."""

    def test_dna_parsing(self):
        """Test parsing of a DNA minimal MEME file."""
        with open(minimal_dna_path) as f:
            motifs_list = read(f, "minimal") 

        self.assertEqual(len(motifs_list), 3) 
        motif = motifs_list[0]
        self.assertEqual(motif.name, "KRP")
        self.assertEqual(motif.length, 19)
        self.assertEqual(motif.consensus[:5], "TGTGA")

    def test_rna_parsing(self):
        """Test parsing of an RNA minimal MEME file."""
        with open(minimal_rna_path) as f:
            motifs_list = read(f, "minimal") 

        self.assertEqual(len(motifs_list), 3)
        motif = motifs_list[0]
        self.assertEqual(motif.name, "KRP") 
        self.assertEqual(motif.length, 19)
        self.assertEqual(motif.consensus[:5], "TGTGA")


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
        self.assertTrue(self.dist.interval > 0.0)
        self.assertTrue(self.dist.ic is not None) 
        self.assertTrue(len(self.dist.mo_density) > 0)
        
    def test_threshold_fpr(self):
        """Test calculation of threshold based on False Positive Rate (FPR)."""
        threshold_low_fpr = self.dist.threshold_fpr(0.01)
        threshold_high_fpr = self.dist.threshold_fpr(0.5)

        self.assertTrue(threshold_low_fpr > threshold_high_fpr)
        self.assertTrue(self.pssm.min < threshold_low_fpr < self.pssm.max)


if __name__ == '__main__':
    # Use unittest.main() to run all tests
    unittest.main(argv=sys.argv[:1], exit=False)
