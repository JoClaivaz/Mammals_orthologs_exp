# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170512, last modification 171008
Extraction of fasta sequence amongst orthologs pairs between specie 1 and specie 2 
(obtained from OMA, current version available at http://omabrowser.org/oma/current/), in function of the pairwise_ortholog file.

Command line
python extract_ortholog_fasta_sequence.py list_considered_species pair_ortho_path sequence_file_input regexp_pair regexp_output
"""

###Function
import os

def extract_ortholog_fasta_sequence (list_considered_specie = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN'], 
                                     pair_ortho_path = 'D:/UNIL/Master/Master_Project/Data/OMA/', 
                                     sequence_file_input = 'D:/UNIL/Master/Master_Project/Data/OMA/fasta_seq/oma-seqs.fa',
                                     regexp_pair = 'pairwise_',
                                     regexp_output = 'protein_sequence_'):
    
    #store all orthologs of considered species into dictionary
    dict_species = {}
    
    for considered_species in list_considered_specie:
        dict_species[considered_species] = {}
    
    pair_ortholog_file = os.listdir(pair_ortho_path)
    
    for PairOrtho in pair_ortholog_file:
        
        if regexp_pair in PairOrtho:
            pair_ortholog = open(pair_ortho_path + PairOrtho, 'r')
            
            for considered_pair in pair_ortholog:
                dict_species[considered_pair.split('\t')[0].rstrip('0123456789')][considered_pair.split('\t')[0]] = ''
                dict_species[considered_pair.split('\t')[1].rstrip('0123456789')][considered_pair.split('\t')[1]] = ''
            
            pair_ortholog.close()
    #
    
    write_switch = False
    sequence_file = open(sequence_file_input, 'r')
    
    #store sequence fasta into dictionnary
    for considered_sequence in sequence_file:
        
        if write_switch and '>' in considered_sequence:
            write_switch = False
            dict_species[protein_tmp.rstrip('0123456789')][protein_tmp] = seq_tmp
        
        elif write_switch:
            seq_tmp += considered_sequence 
        
        try:
            dict_species[considered_sequence.split('> ')[1].replace('\n', '').rstrip('0123456789')][considered_sequence.split('> ')[1].replace('\n', '')] = ''
            write_switch = True
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
        all_ortholog = open('%s%s%s' % (pair_ortho_path, regexp_output, species), 'w')
        
        for sequence in dict_species[species].keys():
            all_ortholog.write('>%s\n%s\n' % (sequence, dict_species[species][sequence]))
        
        all_ortholog.close()
    #
    
#Run the function / WINDOWS
#all species, unique orthologs (1:1) and paralogs
#extract_ortholog_fasta_sequence()
    
#Command line:
#python extract_ortholog_fasta_sequence.py list_considered_species pair_ortho_path sequence_file_input regexp_pair regexp_output 
import sys
extract_ortholog_fasta_sequence(list_considered_species = sys.argv[1],
                                pair_ortho_path = sys.argv[2], 
                                sequence_file_input = sys.argv[3],
                                regexp_pair = sys.argv[4],
                                regexp_output = sys.argv[5])