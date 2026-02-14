# Affine Gapped Global Sequence Alignment2

This repository contains a Python implementation of global DNA sequence alignment using affine gap penalties (Needleman–Wunsch algorithm with gap opening and gap extension scoring).

Files:
- script.py – Python implementation of affine global alignment
- seq1.fa – First input DNA sequence (FASTA format)
- seq2.fa – Second input DNA sequence (FASTA format)
- output_alignment.txt – Output file containing alignment score and optimal alignment

How to Run:
python script.py seq1.fa seq2.fa 2 -1 -5 -1

Parameters:
 2   → Match score
-1  → Mismatch penalty
-5  → Gap opening penalty
-1  → Gap extension penalty

The program reads two FASTA sequences, performs global alignment using affine gap scoring, and writes the optimal alignment and final score to output_alignment.txt.

Time Complexity: O(n × m)

Author: Sanskruti K
BBL-434 Bioinformatics
