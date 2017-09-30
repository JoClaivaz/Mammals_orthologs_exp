#### Mammals_orthologs_exp

###Script folder
##Domain_modification_inference subfolder
*'extract_ortholog_fasta_sequence.py'
extraction of the amino acid fasta sequence in function of the considered ortholog pairs, write files with considered protein sequence.
Input: 'pairwise_ortholog_HUMAN_SPECIE2' and 'oma-seqs.fa'
Output: 'protein_sequence_HUMAN' and 'protein_sequence_SPECIE2'

###Pipeline
##Domain_modification_inference
1. recover fasta sequence from http://omabrowser.org/oma/current/
2. recover homology relationship from http://omabrowser.org/oma/genomePW/ 
3. run python script 'extract_ortolog_fasta_sequence.py'
