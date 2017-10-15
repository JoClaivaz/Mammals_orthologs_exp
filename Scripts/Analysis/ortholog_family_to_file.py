# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 16:40:58 2017

@author: Claivaz
Write one file per ortholog family (only consideration human ortholog)
whole considered ortholog families are present in no_modif or domain_modif group for the different pairwise species comprisons 
"""
import os
path_domain_architecture_files = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/'
list_files = os.listdir(path_domain_architecture_files)

#extract into dictionnary the ortholog pairs present in the different considered domain architecture files
ortholog_species = {}
for file_tmp in list_files:
    if 'ortholog' in file_tmp:
        if 'HUMAN' in file_tmp and 'domain_nomodif' in file_tmp or 'HUMAN' in file_tmp and 'domain_loss' in file_tmp:
            consider_file = open(path_domain_architecture_files + file_tmp, 'r')
            
            for consider_pair in consider_file:
                if 'HUMAN' in consider_pair.split('\t')[0]:
                    try:
                        ortholog_species[consider_pair.split('\t')[0]].append(consider_pair.split('\t')[2])
                    except:
                        ortholog_species[consider_pair.split('\t')[0]] = [consider_pair.split('\t')[2]]
            
                else:
                    try:
                        ortholog_species[consider_pair.split('\t')[2]].append(consider_pair.split('\t')[0])
                    except:
                        ortholog_species[consider_pair.split('\t')[2]] = [consider_pair.split('\t')[0]]            
            
            consider_file.close()
            
#extract into dictionnary DNA sequences
path_sequence_DNA = 'D:/UNIL/Master/Master_Project/Data/OMA/DNA_seq/'
list_files = os.listdir(path_sequence_DNA)
ortholog_DNA = {}

for file_tmp in list_files:
    consider_file = open(path_sequence_DNA + file_tmp, 'r')
    first_tmp = True

    for line_tmp in consider_file:

        if '>' in line_tmp and not first_tmp :
            ortholog_DNA[gene_tmp] = sequence_tmp
            gene_tmp = line_tmp.replace('> ', '').replace('\n', '')
            sequence_tmp = ''

        elif first_tmp:
            first_tmp = False
            gene_tmp = line_tmp.replace('> ', '').replace('\n', '')
            sequence_tmp = ''

        else:    
            sequence_tmp = sequence_tmp + line_tmp
            
    ortholog_DNA[gene_tmp] = sequence_tmp        
    consider_file.close()
    
#write one file per ortholog family
ortholog_number = 1
path_output = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_family/'
for consider_key in ortholog_species.keys():
    family_file = open('%sortholog_family_%i' % (path_output, ortholog_number), 'w')
    ortholog_number += 1
    family_file.write('> %s\n%s\n' % (consider_key, ortholog_DNA[consider_key]))
    
    for consider_ortho in ortholog_species[consider_key]:
        family_file.write('> %s\n%s\n' % (consider_ortho, ortholog_DNA[consider_ortho]))
    
    
    family_file.close()