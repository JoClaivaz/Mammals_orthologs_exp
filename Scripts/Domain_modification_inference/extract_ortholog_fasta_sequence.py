# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170512
Extraction of fasta sequence amongst orthologs pairs between specie 1 and specie 2 
(obtained from OMA, current version available at http://omabrowser.org/oma/current/), in function of the pairwise_ortholog file, only considering 
subset of one parameter (e.g. '1:1', unique ortholog, ':' all pair).
"""

###Function
import os

def extract_ortholog_fasta_sequence (list_considered_specie = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX'], 
                                     pair_ortho_path = 'D:/UNIL/Master/Master_Project/Data/OMA/', 
                                     sequence_file_input = 'D:/UNIL/Master/Master_Project/Data/OMA/fasta_seq/oma-seqs.fa',
                                     regexp_pair_ortho = 'pairwise_ortholog',
                                     regexp_output = 'protein_sequence_',
                                     central_specie = 'HUMAN',
                                     parameter = '1:1'):
    
    #store all orthologs of considered species into dictionary
    dict_species = {}
    dict_species[central_specie] = {} 
    for considered_species in list_considered_specie:
        dict_species[considered_species] = {}
    pair_ortholog_file = os.listdir(pair_ortho_path)
    for PairOrtho in pair_ortholog_file:
        if regexp_pair_ortho in PairOrtho and PairOrtho.split('HUMAN_')[1].replace('.txt', '') in list_considered_specie:
            pair_ortholog = open(pair_ortho_path + PairOrtho, 'r')
            for considered_pair in pair_ortholog:
                if parameter in considered_pair:
                    dict_species[considered_pair.split('\t')[0].rstrip('0123456789')][considered_pair.split('\t')[0]] = ''
                    dict_species[considered_pair.split('\t')[1].rstrip('0123456789')][considered_pair.split('\t')[1]] = ''
            pair_ortholog.close()
    #
    
    write = False
    sequence_file = open(sequence_file_input, 'r')
    
    #store sequence fasta into dictionnary
    for considered_sequence in sequence_file:
        if write and '>' in considered_sequence:
            write = False
            dict_species[protein_tmp.rstrip('0123456789')][protein_tmp] = seq_tmp
        elif write:
            seq_tmp += considered_sequence 
        try:
            dict_species[considered_sequence.split('> ')[1].replace('\n', '').rstrip('0123456789')][considered_sequence.split('> ')[1].replace('\n', '')] = ''
            write = True
            seq_tmp = ''
            protein_tmp = considered_sequence.split('> ')[1].replace('\n', '')
        except KeyError:
            pass
        except IndexError:
            pass 
    dict_species[protein_tmp.rstrip('0123456789')][protein_tmp] = seq_tmp
    sequence_file.close()
    #
    
    #write file
    for species in dict_species.keys():
        unique_ortholog = open('%s%s%s' % (pair_ortho_path, regexp_output, species), 'w')
        for sequence in dict_species[species].keys():
            unique_ortholog.write('>%s\n%s\n' % (sequence, dict_species[species][sequence]))
        unique_ortholog.close()
    #
    
#Run the function / WINDOWS
#all species, unique orthologs (1:1)
extract_ortholog_fasta_sequence()