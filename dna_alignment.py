#!/usr/bin/env python3
"""
dna_alignment.py - A simple tool for pairwise DNA sequence alignment.

This script performs a basic global alignment of two DNA sequences using a
simplified scoring system. It calculates similarity, identifies matches/mismatches,
and provides insights into evolutionary relationships or functional elements.
Suitable for educational purposes.

Author: [Swarnendu Das]
Date: April 04, 2025
"""

def validate_dna_sequence(sequence):
    """
    Validate that a DNA sequence contains only A, T, G, C.

    Args:
        sequence (str): Input DNA sequence

    Raises:
        ValueError: If the sequence contains invalid bases
    """
    valid_bases = set('ATGC')
    sequence_upper = sequence.upper()
    invalid_bases = set(sequence_upper) - valid_bases
    
    if invalid_bases:
        raise ValueError(f"Invalid bases found: {invalid_bases}. Only A, T, G, C are allowed.")

def align_sequences(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    """
    Perform a global alignment of two DNA sequences using a simple scoring system.

    Args:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        match_score (int): Score for a match (default: 1)
        mismatch_penalty (int): Penalty for a mismatch (default: -1)
        gap_penalty (int): Penalty for a gap (default: -2)

    Returns:
        tuple: (aligned_seq1, aligned_seq2, similarity_percentage)
    """
    validate_dna_sequence(seq1)
    validate_dna_sequence(seq2)
    
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    # Initialize the scoring matrix
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # Fill first row and column with gap penalties
    for i in range(rows):
        matrix[i][0] = i * gap_penalty
    for j in range(cols):
        matrix[0][j] = j * gap_penalty
    
    # Fill the matrix
    for i in range(1, rows):
        for j in range(1, cols):
            match = matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = matrix[i-1][j] + gap_penalty
            insert = matrix[i][j-1] + gap_penalty
            matrix[i][j] = max(match, delete, insert)
    
    # Traceback to get aligned sequences
    aligned1, aligned2 = [], []
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + gap_penalty:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    
    # Calculate similarity
    matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != '-')
    total_length = len(aligned1)
    similarity_percentage = (matches / total_length) * 100 if total_length > 0 else 0
    
    return aligned1, aligned2, similarity_percentage

def print_alignment(aligned_seq1, aligned_seq2, similarity):
    """
    Print the aligned sequences with match/mismatch indicators.

    Args:
        aligned_seq1 (str): First aligned sequence
        aligned_seq2 (str): Second aligned sequence
        similarity (float): Similarity percentage
    """
    match_line = ''.join('|' if a == b and a != '-' else ' ' for a, b in zip(aligned_seq1, aligned_seq2))
    print(f"Sequence 1: {aligned_seq1}")
    print(f"            {match_line}")
    print(f"Sequence 2: {aligned_seq2}")
    print(f"Similarity: {similarity:.2f}%")
    print(f"Matches: {sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')}")
    print(f"Mismatches: {sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != '-' and b != '-')}")
    print(f"Gaps: {aligned_seq1.count('-') + aligned_seq2.count('-')}\n")

def main():
    """Main function to demonstrate DNA sequence alignment."""
    # Example sequences
    seq1 = "ATGGCCATAGCTAGCT"
    seq2 = "ATGCGTAGCTAGCTAGC"
    
    try:
        # Align sequences
        aligned_seq1, aligned_seq2, similarity = align_sequences(seq1, seq2)
        
        # Print results
        print("DNA Sequence Alignment Results:")
        print_alignment(aligned_seq1, aligned_seq2, similarity)
        
        # Test with invalid sequence
        print("Testing invalid sequence:")
        align_sequences("ATGCX", seq2)
        
    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
