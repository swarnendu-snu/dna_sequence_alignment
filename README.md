# DNA Sequence Alignment Tool

A simple Python script for pairwise DNA sequence alignment to identify similarities and differences. This tool uses a basic global alignment algorithm with a scoring system, suitable for educational purposes to learn about evolutionary relationships and functional elements.

## Features
- **Pairwise Alignment**: Aligns two DNA sequences globally using a simplified Needleman-Wunsch approach.
- **Similarity Calculation**: Computes the percentage of similarity between sequences.
- **Match/Mismatch Detection**: Identifies matches, mismatches, and gaps in the alignment.
- **Error Handling**: Validates sequences to ensure only A, T, G, C bases are used.

## Limitations
- This is not a replacement for advanced tools like BLAST, which use heuristic methods and large databases.
- Best for short sequences due to computational simplicity.

## Requirements
- Python 3.x (No external libraries required)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/[your-username]/dna_sequence_alignment.git
