#### Mammals_orthologs_exp
This folder is synchronized with our github repo : https://github.com/JoClaivaz/Mammals_homologs_exp.git
This repo is public.
Scripts folder contains all scripts we created for data analysis.

#Aim
Study the effect of domain architecture modification on protein function in mammals evolution.

#Domain Inference
Considered orthologs:
  -Unique orthologs '1:1' from OMA
  -At least one domain inferred by pfam
  -Only one domain modification / no domain modification (control)
  -The modified domain is not repeteted in the protein structure
  -Results cross with domainDIFF (Carsten program, University Mûnster)


7 pairwise species comparisons are considered HUMAN against (Specie2):
  -BOVIN
  -GORGO
  -MACMU
  -MONDO
  -MOUSE
  -PANTR
  -PIGXX
  

###Script folder
##Domain_modification_inference subfolder
*'extract_ortholog_fasta_sequence.py'
extraction of the amino acid fasta sequence in function of the considered ortholog pairs, write files with considered protein sequence.
Input: 'pairwise_ortholog_HUMAN_SPECIE2' and 'oma-seqs.fa'
Output: 'protein_sequence_HUMAN' and 'protein_sequence_SPECIE2'

*'domainDIFF_output_organization_and_execution.py'
create each pairwise file with domain inferred from pfmascan necessary for running domainDiff program (from Carsten, university Münster)
Input: 'Specie2_domain' and 'HUMAN_domain' (output from pfamscan)
Output: 'ortholog_HUMAN_Specie2_domain' and 'ortholog_HUMAN_Specie2_domain_modifications' (from script and domainDiff)

###Pipeline
##Domain_modification_inference
1. recover fasta sequence from http://omabrowser.org/oma/current/
2. recover homology relationship from http://omabrowser.org/oma/genomePW/ 
3. run python script 'extract_ortolog_fasta_sequence.py'
4. run pfamscan perl tool on 'protein_sequence_HUMAN' and 'protein_sequence_SPECIE2' > 'Specie2_domain' and 'HUMAN_domain'
4. run python script 'domainDIFF_output_organization_and_execution.py'
