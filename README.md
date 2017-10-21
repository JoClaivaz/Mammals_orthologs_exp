#### Mammals_orthologs_exp
This folder is synchronized with our github repo : https://github.com/JoClaivaz/Mammals_homologs_exp.git
This repo is public.
Scripts folder contains all scripts we created for data analysis.

#Aim
Study the effect of domain architecture modification on protein function in mammals evolution.
#General HYP
Domain architecture modification has an effect on tissue expression pattern in mammals evolution.

#Domain Inference
Considered orthologs:
  -Unique orthologs '1:1' from OMA
  -At least one domain inferred by pfam
  -Only one domain modification / no domain modification (control)
  -The modified domain is not repeteted in the protein structure
Considered paralogs:
  -Inferred from OMA orthology
  -At least one domain inferred by pfam
  -Only one domain modification / no domain modification (control)
  -The modified domain is not repeteted in the protein structure

9 pairwise species comparisons are considered. Composed of pair of (SPECIE):
	-BOVIN
	-GORGO
	-HUMAN
	-MACMU
	-MONDO
	-MOUSE
	-PANTR
	-PIGXX
	-RATNO


###Script folder
##Domain_modification_inference subfolder
#AIM
Infer pair without domain modification or one domain loss
Filter out domain loss which are repeated in domain architecture, known as less informative and less susceptible to present expresssion pattern shift. 

#'extract_ortholog_pair_OMA.py'
Conserved unique ortholog (1:1) from OMA and infer the different paralog for each species from the not unique ortholog pairs present in OMA. 
A paralog family number is associated in function of which species gene is considered for the paralog inference (WARNING a geneID can have different paralog family; unique paralog family per gene identifier is considered in 'plot_results_tspec.R', can be done before if it accepted)
Input: 'oma-pairs.txt' (recovered from OMA) and list of considered species
Output: 'pairwise_ortholog_SPECIE1_SPECIE2' and 'pairwise_paralog_SPECIE' 
Command line: python extract_ortholog_pair_oma.py list_considered_species pair_ortho_path sequence_file_input regexp_pair regexp_output

#'extract_ortholog_fasta_sequence.py'
extraction of the amino acid fasta sequence in function of the considered ortholog and paralog pairs, store in files considered protein sequences or each species.
Input: 'pairwise_ortholog_SPECIE1_SPECIE2', 'pairwise_paralog_SPECIE' (both called by the combination of 'list_considered_species', 'pair_ortho_path', 'regexp_pair'), 'oma-seqs.fa' (recovered from OMA) and list of considered species
Output: 'protein_sequence_SPECIE'
Command line: python extract_ortholog_fasta_sequence.py list_considered_species pair_ortho_path sequence_file_input regexp_pair regexp_output

#'extract_ortholog_modification_group.py'
extraction of pair which have at least one domain inferred by pfamscan (both/partial/none)
extraction of pair without domain modification (control group in analysis, same domain architecture sequence), complex domain modification (with the same number of domain in both gene of a pair but not the same sequence or more than 1 domain difference), domain modification with 1 domain difference (group of interest, one domain of difference (not consider the kind of domain in this step)). All considered pair have at least one domain inferred by pfamscan.
extraction of pair with one domain modification not involved in domain repetition
Inference of the position of the loss (f-1, b-1, int-1; N-terminal, C-terminal, internal respectively) in the no repeated pair
Input: 'SPECIE_domain' (pfamscan output, called by 'Domain_file_in_path' argument), 'pairwise_ortholog_SPECIE1_SPECIE2' (called by the combination of 'list_mammals' and 'PairOrtho_in_path_with_regexp' arguments)
Output: 'putative_ortholog_SPECIE1_SPECIE2_domain_loss', 'ortholog_SPECIE1_SPECIE2_partial', 'ortholog_SPECIE1_SPECIE2_both', 'ortholog_SPECIE1_SPECIE2_none', 'ortholog_SPECIE1_SPECIE2_domain_nomodif', 'ortholog_SPECIE1_SPECIE2_domain_complex_modif', 'ortholog_SPECIE1_SPECIE2_1_domain_modif', 'ortholog_SPECIE1_SPECIE2_1_domain_repeated' and 'ortholog_SPECIE1_SPECIE2_1_domain_notrepeated'.
N.B.:after each gene in the output file, the domain length is indicated (2nd and 4th columns)
Command line: python extract_ortholog_modification_group.py list_mammals PairOrtho_in_path_with_regexp Domain_file_in_path Domain_file_out_with_regexp

#'extract_paralog_modification_group.py'
extraction of pair which have at least one domain inferred by pfamscan (both/partial/none)
extraction of pair without domain modification (control group in analysis, same domain architecture sequence), complex domain modification (with the same number of domain in both gene of a pair but not the same sequence or more than 1 domain difference), domain modification with 1 domain difference (group of interest, one domain of difference (not consider the kind of domain in this step)). All considered pair have at least one domain inferred by pfamscan.
extraction of pair with one domain modification not involved in domain repetition
Inference of the position of the loss (f-1, b-1, int-1; N-terminal, C-terminal, internal respectively) in the no repeated pair
Input: 'SPECIE_domain' (pfamscan output, called by 'Domain_file_in_path' argument), 'pairwise_paralog_SPECIE' (called by the combination of 'list_mammals' and 'PairPara_in_path_with_regexp' arguments)
Output: 'putative_paralog_SPECIE_domain_loss', 'paralog_SPECIE_partial', 'paralog_SPECIE_both', 'paralog_SPECIE_none', 'paralog_SPECIE_domain_nomodif', 'paralog_SPECIE_domain_complex_modif', 'paralog_SPECIE_1_domain_modif', 'paralog_SPECIE_1_domain_repeated' and 'paralog_SPECIE_1_domain_notrepeated'.
N.B.:after each gene in the output file, the domain length is indicated (2nd and 4th columns) and the 6th column contains the paralog group
Command line: python extract_paralog_modification_group.py list_mammals PairPara_in_path_with_regexp Domain_file_in_path Domain_file_out_with_regexp

#'ortho_domain_modification_analysis.py', 'para_domain_modification_analysis.py' and 'barplot_result_domain_modification.R'
number of gene available analysis in the different files produced through the domain architecture modification inference step, and use R script to obtained plot
Input: all the files produced for 'Domain_modification_inference' for the ortholog part
Ouput: barplot visualization


##Analysis subfolder
#'tissue_availability_mammals_bgee.py' 
extract unique Anatomical entity name and Stage name for each expression dataset from Bgee, allowing the choice of the species
Input: * all expression data from Bgee
Output: 'tissue_available_bgee'

#'extraction_state_expression_file.py'
parsed expression files in function of the considered species and tissues & the results of domain rearrengement (control group / 1 domain lost not repeated in termini part of the proten). 
Formate the outputfile for R analysis (by column): 'ExperimentID\tLibraryID\tGeneID\tAnatomicalEntityName\tStageName\tSex\tDomainStatus\t%s_homolog\tFPKM\n' % (central_species) 
Input: all dataset of 'ortholog_Specie1_Specie2_domain_nomodif', 'putative_ortholog_Specie1_Specie2_domain_loss', 'paralog_SPECIE_domain_nomodif', 'putative_paralog_SPECIE_domain_loss' and expression files from Bgee, 'oma-ensembl' from OMA, list of considered species and tissues
Output: 'SPECIE_expression_parsed'
Command line: python extraction_state_expression_file.py considered_tissues considered_species control_modification_pairs_folder_path  path_expression_dir oma_ensembl_converter output_file_path

#'Tspec_inference.R'
Inference of tissue specificity estimators for ortholog and paralog pairs
Input: 'SPECIE_expression_parsed', 'putative_paralog_SPECIE_domain_loss', 'putative_ortholog_SPECIE_domain_loss', 'paralog_SPECIE_domain_nomodif' and 'ortholog_SPECIE1_SPECIE2_1_domain_modif_domain_nomodif'
Output: store tables for 'Tspec_analysis.R'

#'Tspec_analysis.R'
Study the effect of paralog and ortholog domain modifications on Tspec.
Organization of paralog: reference paralog is defined by the most represented gene for a given paralog group which is whether the longer one or the most expressed one in any tissues (TPM values)
HYP assessed:
	*effect of domain modification on Tspec values for a given species
	*effect of domain modification on correlation between pairwise species comparisons
	*effect of domain modification on Tspec factors shift for a given pairwise species comparison
	*effect of domain architecture length on Tspec values
	*effect of the position where the modification occured
Input: output tables from 'Tspec_inference.R' and 'SPECIE_expression_parsed'

#'plot_results_tspec.R'
Create pdf files containing correaltion results in function of status or kind of modifications in different ortholog or paralog datasets. 
Consideration of unique paralog pair, having an unique Paralog Group identifier: longest and/or most expresssed gene of pair in any considered tissues.
Plot:
	*tspec Specie1 function of tspec Specie2
	*tspec Specie1 function of tspec Specie2 and function of the modification group
Input: output tables from 'Tspec_inference.R' and 'SPECIE_expression_parsed'
Output: pdf contiainig plots

#'cor_dif_test.R'
Function allowing the process of exact fisher test on the correlation amongst specie 1 and specie2, and other one allows to get the number of tissues considered for a given analysis
Input: output tables from 'Tspec_inference.R' and 'SPECIE_expression_parsed'

******NOT FINISH********
#'extract_DNA_sequence.py'
create one file containing all the cDNA sequence specific to one specie extracted from 'eukaryotes.cdna.fa'
Input: 'eukaryotes.cdna.fa' and list of considered species
Output: 'SPECIE_dna_seq'

#'ortholog_family_to_file.py' and 'paralog_family_to_file.py'
create one file per ortholog or paralog family, contining cDNA sequences. Considered pairs are present in control or domain modification considered groups.
Input: 'ortholog_Specie1_Specie2_domain_nomodif', 'paralog_SPECIE_domain_nomodif', 'putative_ortholog_SPECIE1_SPECIE2_domain_loss', 'putative_paralog_SPECIE_domain_loss' and 'SPECIE_dna_seq'
Output: 'ortholog_family_#' and 'SPECIE_family_#' (paralog files)

#'kaks_analysis.R'
compute the selective pressure amongst homolog genes based on the calculation of kaks
Input: 'ortholog_family_#' and 'SPECIE_family_#'
********************


###Pipeline
##Domain_modification_inference
	1. recover protein fasta sequence ('oma-seqs.fa') and pairwise orthology ('oma-pairs.txt', read with 'zcat') from http://omabrowser.org/oma/current/
	2. run python script 'extract_ortholog_pair_OMA.py'
	3. run python script 'extract_ortolog_fasta_sequence.py'
	4. run pfamscan perl tool on 'protein_sequence_SPECIE' > 'Specie_domain'
	5. run 'extract_ortholog_modification_group.py'
	6. Use 'ortho_domain_modification_analysis.py' and 'barplot_result_domain_modification.R'

##Analysis
	1. recover expression data from http://bgee.org/
	2. recover a table of conversion between OMA and Bgee identifier ('oma-ensembl') from http://omabrowser.org/oma/current/
	3. run 'tissue_availability_mammals_bgee.py' (not mendatory if the species are selected and the states are known)
	4. run 'extraction_state_expression_file.py'
	5. run 'Tspec_inference.R'
	6. use 'Tspec_analysis.R', 'plot_results_tspec.R' and 'cor_dif_test.R'
	7.
	*****NOT FINISH****
	8. recover cDNA fasta sequences ('eukaryotes.cdna.fa') from http://omabrowser.org/oma/current/
	9. run 'extract_DNA_sequence.py'
	10. run 'ortholog_family_to_file.py' and 'paralog_family_to_file.py'
	11. use 'kaks_analysis.R'