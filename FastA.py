import regex as re
import numpy as np
#read the fasta file with the sequences of all yeast genes and their names
def readmultifasta(file):
    f = open(file, 'r')
    sequences = []
    genes = []
    for line in f:
        x=re.match(">", line)
        if x == None:
            sequences.append(line.rstrip())
        else:
            genes.append(line.rstrip())
            
    return  genes, sequences

#read the fasta file with the query sequence 
def readqueryfile(file):
    sequence = ""
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence

#dictionary with keys the kmers of the sequence and items their indexes in the sequence
def Kmers_and_indexes(sequence, k):
    kmers_and_indexes = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer not in kmers_and_indexes:
            kmers_and_indexes[kmer] = [i]
        else:
            kmers_and_indexes[kmer] += [i] 
    return kmers_and_indexes

def diagonal_dictionary(sequence_indexes_dictionary,querykmers_indexes):
    diag_dict = {}
    for kmer in querykmers_indexes.keys():
        if kmer in sequence_indexes_dictionary.keys():
            i_s = querykmers_indexes[kmer]
            p_s = sequence_indexes_dictionary[kmer]
            for j in range(len(i_s)):
                for k in range(len(p_s)):
                    d = i_s[j] - p_s[k]
                    if d in diag_dict.keys():
                        diag_dict[d] += 1
                    else:
                        diag_dict[d] = 1
    return diag_dict

seqs=readmultifasta("all_yeast_genes_minplus1k.fa")[1]
genenames = readmultifasta("all_yeast_genes_minplus1k.fa")[0]
query = readqueryfile("query.fa")

querykmers_indexes = Kmers_and_indexes(query,5)
Matching_Sequences = []
Scores_of_Matching_Sequences = []
minimum_score = len(query) // 2
joined_diagonals_dict = {}
for num in range(1,len(genenames)):
    seqkmers_indexes = Kmers_and_indexes(seqs[num],5)
    diag_dict = diagonal_dictionary(seqkmers_indexes, querykmers_indexes)
    matched_diag_dict = {d: s for d, s in diag_dict.items() if s >= 20}
    matched_diag_dict = dict(sorted(matched_diag_dict.items(), key=lambda x: x[0]))
    current_diagonal = []
    for d, count in matched_diag_dict.items():
        if not current_diagonal:
            current_diagonal.append(d)
            score = count + 4 
        elif d - current_diagonal[-1] <= 3:  # gap <= 3
            current_diagonal.append(d)
            score += count + 4  
        else:
            if score  >= minimum_score:
                joined_diagonals_dict[current_diagonal] = score
            current_diagonal = [d]
    if len(current_diagonal) ==1 : #check the last element (d) of the current_diagonal that is not checked
        if matched_diag_dict[current_diagonal[0]]+4 >= minimum_score:
            joined_diagonals_dict[current_diagonal[0]] = matched_diag_dict[current_diagonal[0]]+4
    else:
        for d in current_diagonal:  #In case g < 3 and the score hasn't been checked because there hasn't been a gap over 3 at the end in order for the score to be checked
            if score >= minimum_score:  
                joined_diagonals_dict[current_diagonal] = score
    if joined_diagonals_dict:
        Matching_Sequences.append(num)
        Scores_of_Matching_Sequences.append(max(joined_diagonals_dict.values()))
        joined_diagonals_dict = {}
        
index = Scores_of_Matching_Sequences.index(max(Scores_of_Matching_Sequences))
print(f"Highest matching sequence: {Matching_Sequences[index]}, gene name: {genenames[Matching_Sequences[index]]} , matching score: {Scores_of_Matching_Sequences[index]}")
