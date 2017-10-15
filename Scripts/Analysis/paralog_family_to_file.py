# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 17:38:38 2017

@author: Claivaz

Write one file per paralog family
whole considered paralog families are present in no_modif or domain_modif group
"""
def paralog_family_to_file(considered_specie,
                           path_domain_architecture_files = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/',
                           path_sequence_DNA = 'D:/UNIL/Master/Master_Project/Data/OMA/DNA_seq/',
                           path_output = 'D:/UNIL/Master/Master_Project/Data/OMA/paralog_family/'):

    import os
    list_files = os.listdir(path_domain_architecture_files)
    
    #extract into dictionnary the paralog pairs present in the different considered domain architecture files
    paralog_species = {}
    for file_tmp in list_files:
        if 'paralog' in file_tmp and considered_specie in file_tmp:
            if 'domain_nomodif' in file_tmp or 'domain_loss' in file_tmp:
                consider_file = open(path_domain_architecture_files + file_tmp, 'r')
                
                for consider_pair in consider_file:
                    try:
                        if consider_pair.split('\t')[0] not in paralog_species[consider_pair.split('\t')[5]]:
                            paralog_species[consider_pair.split('\t')[5]].append(consider_pair.split('\t')[0])
                            
                        if consider_pair.split('\t')[2] not in paralog_species[consider_pair.split('\t')[5]]:
                            paralog_species[consider_pair.split('\t')[5]].append(consider_pair.split('\t')[2])
                    
                    except:
                        paralog_species[consider_pair.split('\t')[5]] = [consider_pair.split('\t')[0]]
                        paralog_species[consider_pair.split('\t')[5]].append(consider_pair.split('\t')[2])
                        
                consider_file.close()
                
    #extract into dictionnary DNA sequences
    paralog_DNA = {}
    consider_file = open('%s%s_dna_seq' % (path_sequence_DNA, considered_specie), 'r')
    first_tmp = True

    for line_tmp in consider_file:

        if '>' in line_tmp and not first_tmp :
            paralog_DNA[gene_tmp] = sequence_tmp
            gene_tmp = line_tmp.replace('> ', '').replace('\n', '')
            sequence_tmp = ''

        elif first_tmp:
            first_tmp = False
            gene_tmp = line_tmp.replace('> ', '').replace('\n', '')
            sequence_tmp = ''

        else:    
            sequence_tmp = sequence_tmp + line_tmp
            
    paralog_DNA[gene_tmp] = sequence_tmp        
    consider_file.close()
        
    #write one file per paralog family
    for consider_key in paralog_species.keys():
        family_file = open('%s%s_family_%s' % (path_output, considered_specie, consider_key), 'w')
        
        for consider_paralog in paralog_species[consider_key]:
            family_file.write('> %s\n%s\n' % (consider_paralog, paralog_DNA[consider_paralog]))
        
        
        family_file.close()

#Run function in different considered species
for specie_tmp in ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN']:
    paralog_family_to_file(considered_specie = specie_tmp)