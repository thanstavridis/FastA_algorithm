# FastA_algorithm
This project is a simplified implementation of the FastA algorithm, designed to compare a query DNA sequence against the complete set of coding (gene) sequences of Saccharomyces cerevisiae (baker's yeast). The purpose is to identify significant matches based on word-level similarity using a lightweight scoring system and filtering criteria.

# Objective
To simulate a basic version of the FastA sequence alignment algorithm with the following parameters:

- Scoring Scheme:
+1 for each exact word match (5-mer), 0 for mismatches

- Word size/ k-mer size (k): 5

- Minimum matched diagonal size (m): 20

- Maximum gap allowed when joining regions (g): 3

- Threshold for total score: 50% of the query sequence length

- Outputs:

  - Identification of the best matching sequence based on the matching highest score

  - Identification of the best matching gene based on the highest score

  - Highest matching score

# Data format
- Query sequence: single-sequence FASTA file

- Gene sequences (S. cerevisiae): multi- FASTA file
