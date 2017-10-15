# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 15:49:54 2017

@author: Claivaz
The script extract the considered mRNA sequences for ka/ks analysis
"""

def extract_DNA_sequence_OMA(list_species = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 
                                             'PANTR', 'PIGXX', 'RATNO', 'HUMAN'],
                             dna_sequence_oma_input = 'E:/OMA/eukaryotes.cdna.fa',
                             folder_output = 'D:/UNIL/Master/Master_Project/Data/OMA/DNA_seq/'):
    
    dna_file = open(dna_sequence_oma_input, 'r')
    considered_species = ''
    write_sp = False
    
    for dna_line in dna_file:
        if write_sp and '>' in dna_line:
            if dna_line.replace('> ', '').replace('\n', '').rstrip('0123456789') == considered_species:
                output_file.write(dna_line)
        
            else: 
                write_sp = False
                output_file.close()
            
        elif write_sp:
            output_file.write(dna_line)
        
        if not write_sp and dna_line.replace('> ', '').replace('\n', '').rstrip('0123456789') in list_species:
                considered_species = dna_line.replace('> ', '').replace('\n', '').rstrip('0123456789')
                output_file = open('%s%s_dna_seq' % (folder_output, considered_species), 'w')
                output_file.write(dna_line)
                write_sp = True

    dna_file.close()

#Run function
extract_DNA_sequence_OMA()