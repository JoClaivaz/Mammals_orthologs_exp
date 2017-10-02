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


8 pairwise species comparisons are considered HUMAN against (Specie2):
	-BOVIN
	-GORGO
	-MACMU
	-MONDO
	-MOUSE
	-PANTR
	-PIGXX
	-RATNO


###Script folder
##Domain_modification_inference subfolder
#AIM
Filter out domain loss which are repeated in domain architecture, which are less informative and less susceptible to present expresssion pattern shift
#'extract_ortholog_fasta_sequence.py'
extraction of the amino acid fasta sequence in function of the considered ortholog pairs, write files with considered protein sequence.
Input: 'pairwise_ortholog_HUMAN_SPECIE2.txt' and 'oma-seqs.fa'
Output: 'protein_sequence_HUMAN' and 'protein_sequence_SPECIE2'

#'domainDIFF_output_organization_and_execution.py'
create each pairwise file with domain inferred from pfmascan necessary for running domainDiff program (from Carsten, university Münster)
Input: 'Specie2_domain' and 'HUMAN_domain' (output from pfamscan)
Output: 'ortholog_HUMAN_Specie2_domain' and 'ortholog_HUMAN_Specie2_domain_modifications' (from script and domainDiff)

#'extract_ortholog_modification_group.py'
extraction of pair without domain modification (control group in analysis, same domain architecture sequence), complex domain modification (with the same number of domain in both gene of a pair but not the same sequence or more than 1 domain difference), domain modification with 1 domain difference (group of interest, one domain of difference (not consider the kind of domain in this step)). All considered pair have at least one domain inferred by pfamscan.
Input: 'ortholog_HUMAN_Specie2_domain' and 'pairwise_ortholog_HUMAN_Specie2.txt'
Output: 'ortholog_HUMAN_Specie2_domain_nomodif', 'ortholog_HUMAN_Specie2_domain_complex_modif' and 'ortholog_HUMAN_Specie2_1_domain_modif'.

#'extract_ortholog_repeated_domain.py'
extration of pair with 1 domain modification and not involved in domain repetition
Input: 'ortholog_HUMAN_Specie2_1_domain_modif' and 'ortholog_HUMAN_Specie2_domain' (output pfamscan)
Output: 'ortholog_HUMAN_Species2_1_domain_notrepeated', 'ortholog_HUMAN_Species2_1_domain_repeated', 'ortholog_HUMAN_Species2_1_domain_complex_modif'

#'extract_ortholog_pfamscan_status.py'
extraction of pair which have at least one domain inferred by pfamscan (both/partial/none), usefull for analysis of this step
Input: 'ortholog_HUMAN_Specie2_domain' and 'pairwise_ortholog_HUMAN_Specie2.txt'
Output: 'ortholog_HUMAN_Specie2_both', 'ortholog_HUMAN_Specie2_partial' and 'ortholog_HUMAN_Specie2_none'.

#'extract_ortholog_domain_loss_domainDiff.py'
intersect of ortholog pair with one domain modification in N- or C- termini of the protein inferred by domainDIFF and the loss domain is not part of a repetition inferred by 'extract_ortholog_repeated_domain.py'
Input: 'ortholog_HUMAN_Species2_1_domain_notrepeated' and 'ortholog_HUMAN_Specie2_domain_modifications'
Output: 'final_pair_HUMAN_Specie2_domain_loss'

#'domain_modification_analysis.py' and 'barplot_result_domain_modification.R'
analysis number of gene available in the different files produced through this step, and use R script to obtained plot
Input: all the files produced for 'Domain_modification_inference'
Ouput: barplot visualization

###Pipeline
##Domain_modification_inference
	1. recover fasta sequence from http://omabrowser.org/oma/current/
	2. recover homology relationship from http://omabrowser.org/oma/genomePW/ 
	3. run python script 'extract_ortolog_fasta_sequence.py'
	4. run pfamscan perl tool on 'protein_sequence_HUMAN' and 'protein_sequence_SPECIE2' > 'Specie2_domain' and 'HUMAN_domain'
	5. run python script 'domainDIFF_output_organization_and_execution.py'
	6. run 'extract_ortholog_modification_group.py'
	7. run 'extract_ortholog_repeated_domain.py'
	8. run 'extract_ortholog_pfamscan_status.py'
	9. run 'extract_ortholog_domain_loss_domainDiff.py'
	10. Use 'domain_modification_analysis.py' and 'barplot_result_domain_modification.R'