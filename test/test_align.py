# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    # using instructions above
    aligner = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1) 
    aligner.align(seq1, seq2) # perform alignment to fill in matrices

    # test that matrices are filled with proper shape
    n_seqA = len(seq1)
    n_seqB = len(seq2)
    assert aligner._align_matrix.shape == (n_seqA + 1, n_seqB + 1), f"Alignment matrix shape is incorrect: expected {(n_seqA + 1, n_seqB + 1)}, got {aligner._align_matrix.shape}"
    assert aligner._gapA_matrix.shape == (n_seqA + 1, n_seqB + 1), f"GapA matrix shape is incorrect: expected {(n_seqA + 1, n_seqB + 1)}, got {aligner._gapA_matrix.shape}"
    assert aligner._gapB_matrix.shape == (n_seqA + 1, n_seqB + 1), f"GapB matrix shape is incorrect: expected {(n_seqA + 1, n_seqB + 1)}, got {aligner._gapB_matrix.shape}"
    
    # test that matrices are not all -inf (they should have been filled)
    assert not np.all(np.isinf(aligner._align_matrix)), "Alignment matrix is not filled correctly, all values are -inf" 
    assert not np.all(np.isinf(aligner._gapA_matrix)), "GapA matrix is not filled correctly, all values are -inf"
    assert not np.all(np.isinf(aligner._gapB_matrix)), "GapB matrix is not filled correctly, all values are -inf"

    # test that base cases are initialized correctly
    assert aligner._align_matrix[0, 0] == 0, f"Base case align_matrix[0,0] should be 0, got {aligner._align_matrix[0, 0]}"
    assert aligner._gapA_matrix[0, 0] == -np.inf, f"Base case gapA_matrix[0,0] should be -inf, got {aligner._gapA_matrix[0, 0]}"
    assert aligner._gapB_matrix[0, 0] == -np.inf, f"Base case gapB_matrix[0,0] should be -inf, got {aligner._gapB_matrix[0, 0]}"
    
    # test affine gap penalties for first row and column
    # first row: gap_open + (j-1) * gap_extend
    for j in range(1, len(seq2) + 1):
        expected = -10 + (j - 1) * (-1)
        assert aligner._gapB_matrix[0, j] == expected, f"GapB matrix first row not initialized correctly at column {j}: expected {expected}, got {aligner._gapB_matrix[0, j]}"
    
    # first column: gap_open + (i-1) * gap_extend
    for i in range(1, len(seq1) + 1):
        expected = -10 + (i - 1) * (-1)
        assert aligner._gapA_matrix[i, 0] == expected, f"GapA matrix first column not initialized correctly at row {i}: expected {expected}, got {aligner._gapA_matrix[i, 0]}"
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    aligner = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1) # create aligner
    score, seqA_align, seqB_align = aligner.align(seq3, seq4) # perform alignment and get backtrace results
    
    # test alignment score
    assert score == 18, f"Expected score 18, got {score}"
    
    # test alignment strings
    assert seqA_align == "MAVHQLIRRP", f"Expected seqA 'MAVHQLIRRP', got '{seqA_align}'"
    assert seqB_align == "M---QLIRHP", f"Expected seqB 'M---QLIRHP', got '{seqB_align}'"
    
    # test that alignments have the same length
    assert len(seqA_align) == len(seqB_align), "Aligned sequences should have the same length"




