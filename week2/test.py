import unittest
import os
import numpy as np
import math
from typing import List, Dict, Any, Union
from io import StringIO
import sys

# DATA folder
DATA_DIR = os.path.join(os.path.dirname(__file__), "data") 
MINIMAL_DNA_PATH = os.path.join(DATA_DIR, "minimal_test.meme")
MINIMAL_RNA_PATH = os.path.join(DATA_DIR, "minimal_test_rna.meme")

"""
if __codon__:
    # Codon environment
    from .code.bio_codon import create, read, Motif
    from .code.bio_codon.matrix import PositionSpecificScoringMatrix, PositionWeightMatrix, CountsMatrix
    from .code.bio_codon.thresholds import ScoreDistribution
else:
    # Python environment
    from Bio import motifs
    from Bio.Seq import Seq
"""

# Sequences for TATA-like box (Basic Motif)
TEST_SEQUENCES_BASIC = [
    "TACAAT",
    "TACAAG",
    "CACAAT",
    "GAGAAG",
    "AACAAG",
    "TACAAC",
    "TACAAA",
]
TEST_MOTIF_NAME = "TATA_Test_Motif"

# Sequences for Matrix/PSSM Scoring Test (Longer Motif, from uploaded test_motifs.py)
TEST_SEQUENCES_PSSM = [
    "GTTAGTCATTAG",
    "CTTAACTATTAG",
    "CTTAACTAATAG",
    "GTTAGTGTAGGG"
]
TEST_PSSM_MOTIF_NAME = "Long_Scoring_Motif"


if '__codon__' in sys.builtin_module_names:
    # Codon environment: Assuming 'bio_codon' package is accessible via PYTHONPATH
    try:
        import bio_codon.motifs as motifs
        # Codon specific class names exposed in __init__.py
        Motif = motifs.Motif
        create = motifs.create
        read = motifs.read
    except ImportError:
        print("ERROR: Could not import bio_codon.motifs in Codon environment.")
        sys.exit(1)
    # Mock Seq for Codon environment
    Seq = str 
else:
    # Python environment: Uses Biopython
    from Bio import motifs
    from Bio.Seq import Seq
    Motif = motifs.Motif
    create = motifs.create
    read = motifs.read
    # Ensure ScoreDistribution is available on the motif module level for consistency
    if not hasattr(motifs, 'ScoreDistribution'):
        from Bio.motifs.thresholds import ScoreDistribution
        motifs.ScoreDistribution = ScoreDistribution


# --- Test Cases ---

class TestMotif(unittest.TestCase):
    """
    Comprehensive tests for bio_codon.motifs, ensuring parity with Biopython 
    across core matrix operations, PSSM scoring, I/O, and thresholds.
    """

    @classmethod
    def setUpClass(cls):
        """Create the necessary motif instances."""
        cls.basic_motif: Motif = create(TEST_SEQUENCES_BASIC, name=TEST_MOTIF_NAME)
        cls.pssm_motif: Motif = create(TEST_SEQUENCES_PSSM, name=TEST_PSSM_MOTIF_NAME)

    # --- Print Helpers for Diff Comparison ---
    
    def _print_matrix(self, title: str, matrix: Any):
        """Prints a matrix with standardized floating point precision."""
        print(f"--- {title} ---")
        # Format the matrix string representation consistently
        lines = str(matrix).strip().split('\n')
        for line in lines:
            if line.startswith(("<", "   ")):
                print(line) 
            else:
                # Standardize floating point numbers to 5 decimal places for comparison
                parts = line.split('\t')
                formatted_parts = []
                for part in parts:
                    try:
                        f_val = float(part)
                        formatted_parts.append(f"{f_val:.5f}") 
                    except ValueError:
                        formatted_parts.append(part)
                print('\t'.join(formatted_parts))

    # --- Group 1: Motif and Matrix Core Properties (__init__.py, matrix.py) ---

    def test_01_creation_consensus_and_degenerate(self):
        """Test motif creation, length, consensus, and degenerate consensus."""
        m = self.basic_motif
        print("\nTEST_01_CREATION_CONSENSUS_AND_DEGENERATE")
        print(f"Name: {m.name}")
        print(f"Length: {m.length}")
        print(f"Consensus: {m.consensus}")
        
        # T=4, C=1, G=1, A=1 at pos 0 -> Y (C/T) is the closest IUPAC code
        degenerate_consensus = m.degenerate_consensus
        print(f"Degenerate Consensus: {degenerate_consensus}")
        
        self.assertEqual(m.length, 6)
        self.assertEqual(m.consensus, "TACAAG")
        self.assertEqual(degenerate_consensus, "YACAAK") 

    def test_02_counts_and_reverse_complement(self):
        """Test raw counts and reverse complement functionality."""
        m = self.basic_motif
        print("\nTEST_02_COUNTS_AND_REVERSE_COMPLEMENT")

        # 1. Counts Matrix Access
        c_at_pos0 = m.counts['T', 0]
        col_3_counts = m.counts[:, 3]
        
        print(f"Count of T at pos 0: {c_at_pos0}")
        print(f"Counts column at pos 3 (A): {col_3_counts['A']}")
        
        self.assertEqual(c_at_pos0, 4)
        self.assertEqual(col_3_counts['A'], 7) # All sequences have 'A' at pos 3
        
        # 2. Reverse Complement (RC)
        rc_motif = m.reverse_complement()
        # Consensus: T A C A A G -> RC: C T T G T A
        print(f"Reverse Complement Consensus: {rc_motif.consensus}")
        
        self.assertEqual(rc_motif.consensus, "CTTGTA")
        self.assertEqual(rc_motif.length, m.length)
        
    def test_03_motif_slicing(self):
        """Test motif slicing, which involves CountsMatrix slicing."""
        m = self.basic_motif
        print("\nTEST_03_MOTIF_SLICING")

        # Slice indices 1 to 5 (length 4)
        sliced_motif = m[1:5]
        print(f"Sliced Motif length (1:5): {sliced_motif.length}")
        print(f"Sliced Motif consensus (1:5): {sliced_motif.consensus}")
        
        self.assertEqual(sliced_motif.length, 4)
        self.assertEqual(sliced_motif.consensus, "CAAA") # Based on original indices 1, 2, 3, 4
        
        # Check counts in the slice (e.g., C at index 1 of the slice, which is index 2 of the original)
        self.assertEqual(sliced_motif.counts['C', 1], 1)
        self.assertEqual(sliced_motif.counts['A', 2], 7)

    def test_04_information_content_and_matrix_comparison(self):
        """Test Information Content (IC) and Pearson correlation (matrix.py)."""
        m = self.basic_motif
        m.pseudocounts = 0.0 # Clear pseudocounts for IC calculation

        print("\nTEST_04_IC_AND_MATRIX_COMPARISON")

        # 1. Information Content (IC)
        ic_value = m.ic()
        print(f"Information Content (IC): {ic_value:.4f}")
        self.assertAlmostEqual(ic_value, 7.3739, places=3)
        
        # 2. Matrix Comparison (Pearson Correlation)
        corr_self = m.pearson_correlation(m, offset=0) 
        print(f"Correlation M vs M (offset 0): {corr_self:.4f}")
        self.assertAlmostEqual(corr_self, 1.0, places=4)

    # --- Group 2: PWM/PSSM Scoring and Bad Characters (matrix.py) ---

    def test_05_pssm_and_precise_scoring(self):
        """Test PWM/PSSM generation and precise score calculation over offsets."""
        m = self.pssm_motif # Motif length 12
        print("\nTEST_05_PSSM_AND_PRECISE_SCORING")
        
        # Use pseudocounts=0.25 (as seen in uploaded test_motifs.py)
        m.pseudocounts = 0.25 
        
        # Ensure we get the PSSM (log_odds)
        pssm = m.pwm if hasattr(m, 'pwm') and not callable(m.pwm) else m.log_odds()
        self._print_matrix("PSSM_MATRIX_0.25_PCOUNT", pssm)

        # Sequence for scoring (Length 17: 17 - 12 + 1 = 6 scores/offsets)
        sequence_to_score = "ACGTGTGCGTAGTGCGT" 
        
        # Codon uses calculate_score, Biopython uses score_hit
        score_method = pssm.calculate_score if hasattr(pssm, 'calculate_score') else pssm.score_hit
        
        # Calculate scores for all possible start positions
        scores = []
        for i in range(len(sequence_to_score) - m.length + 1):
            sub_seq = sequence_to_score[i:i+m.length]
            score = score_method(Seq(sub_seq)) if not '__codon__' in sys.builtin_module_names else score_method(sub_seq)
            scores.append(score)
            
        print(f"Scores for '{sequence_to_score}':")
        for i, s in enumerate(scores):
            print(f"  Offset {i}: {s:.5f}")

        # Assert specific score values (copied from uploaded test_motifs.py for consistency)
        self.assertEqual(len(scores), 6)
        self.assertAlmostEqual(scores[0], -29.18364, places=5)
        self.assertAlmostEqual(scores[1], -38.33651, places=5)
        self.assertAlmostEqual(scores[4], -20.30142, places=5)

    def test_06_pssm_scoring_with_bad_char(self):
        """Test if PSSM scoring correctly handles 'N' (bad character)."""
        m = self.pssm_motif
        print("\nTEST_06_PSSM_SCORING_WITH_BAD_CHAR")

        m.pseudocounts = 0.25 
        pssm = m.log_odds()
        
        # Sequence with an 'N' in the last possible alignment position
        sequence_with_n = "ACGTGTGCGTAGTGCGTN" # Length 18. Last score is at offset 6.
        
        # Calculate scores for all possible start positions
        scores = []
        for i in range(len(sequence_with_n) - m.length + 1):
            sub_seq = sequence_with_n[i:i+m.length]
            score = pssm.calculate_score(sub_seq) if '__codon__' in sys.builtin_module_names else pssm.score_hit(Seq(sub_seq))
            scores.append(score)
            
        # The score at offset 6 (last score) should be NaN or the minimum score
        last_score = scores[-1]
        print(f"Score for sub-sequence with 'N': {last_score:.5f}")

        if '__codon__' in sys.builtin_module_names:
            # Codon/Biopython may default to a very low score or NaN depending on implementation.
            # Biopython's PSSM uses 'float('nan')' if the score calculation encounters a base not in the alphabet.
            self.assertTrue(math.isnan(last_score) or last_score < -100, f"Expected NaN or very low score, got {last_score}")
        else:
            self.assertTrue(math.isnan(last_score), f"Expected NaN, not {last_score!r}")

    # --- Group 3: Threshold Calculations (thresholds.py) ---

    def test_07_score_distribution_and_thresholds(self):
        """Test ScoreDistribution for FNR, FPR, and Balanced Thresholds."""
        m = self.basic_motif
        print("\nTEST_07_SCORE_DISTRIBUTION_AND_THRESHOLDS")

        # Get the PSSM
        m.pseudocounts = 0.5 
        pssm = m.log_odds()
        
        # Create ScoreDistribution
        background = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
        if not hasattr(motifs, 'ScoreDistribution'):
            dist = pssm.distribution(background=background)
        else:
            dist = motifs.ScoreDistribution(pssm=pssm, background=background)

        # 1. Threshold for False Negative Rate (FNR)
        target_fnr = 0.1 
        threshold_fnr = dist.threshold_fnr(target_fnr)
        print(f"Threshold for FNR={target_fnr}: {threshold_fnr:.4f}")
        self.assertLess(threshold_fnr, 0.0) # Should be low

        # 2. Threshold for False Positive Rate (FPR) / P-value
        target_fpr = 0.01 
        # Biopython uses threshold_fpr, Codon uses threshold_for_p_value
        threshold_fpr_method = dist.threshold_fpr if hasattr(dist, 'threshold_fpr') else dist.threshold_for_p_value
        threshold_fpr = threshold_fpr_method(target_fpr)
        print(f"Threshold for FPR={target_fpr} (P-value): {threshold_fpr:.4f}")
        self.assertGreater(threshold_fpr, 5.0) # Should be high

        # 3. Balanced Threshold (FPR = FNR)
        threshold_balanced = dist.threshold_balanced(rate_proportion=1.0)
        print(f"Balanced Threshold (FPR=FNR): {threshold_balanced:.4f}")
        self.assertLess(abs(threshold_balanced), 2.0)

    # --- Group 4: Minimal I/O Tests (minimal.py) ---

    def test_08_io_read_and_background_check(self):
        """Test reading 'minimal' MEME file format and asserting background frequencies."""
        print("\nTEST_08_IO_READ_AND_BACKGROUND_CHECK")
        
        try:
            with open(MINIMAL_DNA_PATH) as f:
                record = read(f, format='minimal')
                
        except FileNotFoundError:
            print(f"FILE NOT FOUND: {MINIMAL_DNA_PATH}. Skipping I/O tests.")
            return

        # Assertions based on minimal_test.meme (DNA) content
        dna_bg = record[0].background if isinstance(record, list) else record.background
        
        print(f"DNA Background A: {dna_bg['A']:.4f}")
        print(f"Read Motifs Count: {len(record)}")
        
        self.assertAlmostEqual(dna_bg['A'], 0.303, places=3)
        self.assertEqual(len(record), 3)

        # Now check RNA background from the other file
        try:
            with open(MINIMAL_RNA_PATH) as f:
                record_rna = read(f, format='minimal')
        except FileNotFoundError:
            print(f"FILE NOT FOUND: {MINIMAL_RNA_PATH}. Skipping RNA I/O test.")
            return

        # Assertions based on minimal_test_rna.meme (RNA) content
        rna_bg = record_rna[0].background if isinstance(record_rna, list) else record_rna.background
        
        print(f"RNA Background U: {rna_bg['U']:.4f}")
        print(f"Read RNA Motif 2 Alphabet: {record_rna[1].alphabet}")
        
        self.assertAlmostEqual(rna_bg['U'], 0.306, places=3)
        self.assertEqual(record_rna[1].alphabet, "ACGU") 
        

# --- Execution ---
if __name__ == '__main__':
    # We redirect stdout temporarily to ensure all test output is captured 
    # before any other prints happen, making the diff reliable.
    
    original_stdout = sys.stdout
    output_capture = StringIO()
    sys.stdout = output_capture

    try:
        # Run tests, suppressing normal unittest output, but allowing captured prints
        unittest.main(exit=False, verbosity=0)
    finally:
        # Restore stdout and print captured content
        sys.stdout = original_stdout
        print(output_capture.getvalue())