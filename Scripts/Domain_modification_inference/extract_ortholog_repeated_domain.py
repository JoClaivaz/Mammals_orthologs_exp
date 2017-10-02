# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 14:10:58 2017

@author: Claivaz
Extraction of orthologs pairs between specie 1 and specie 2 having 1 domain modification 
which have 1 domain of difference not repeated, repeated or complex modification. 
"""

def extract_ortholog_repeated_domain(PairOrtho_1domain_in, 
                                     Domain_file_in, 
                                     Domain_file_out,
                                     central_species,
                                     considered_species):
    
    PairOrtho = open(PairOrtho_1domain_in, 'r')
    
    #create gene dictionnary
    spec_domain = {}
    spec_domain[central_species] = {}
    spec_domain[considered_species] = {}
    
    for pair_ortho in PairOrtho:
        spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]] = []
        spec_domain[pair_ortho.split('\t')[2].rstrip('0123456789')][pair_ortho.split('\t')[2]] = []

    PairOrtho.close()
    #
    
    #store domain into dictionnary
    Domain_file = open(Domain_file_in, 'r')
    for ortho_domain_line in Domain_file:
        if ortho_domain_line.replace(' ','\t').split('\t')[0] in spec_domain[central_species].keys() or ortho_domain_line.replace(' ','\t').split('\t')[0] in spec_domain[considered_species].keys():
            ortho_domain_line = ortho_domain_line.replace(' ', '\t')
            while '\t\t' in ortho_domain_line:
                ortho_domain_line = ortho_domain_line.replace('\t\t', '\t')
            spec_domain[ortho_domain_line.split('\t')[0].rstrip('0123456789')][ortho_domain_line.split('\t')[0]].append(ortho_domain_line.split('\t')[6])
     #   
        
    #open file for different group modification    
    Domain_file_notrepeated = open(Domain_file_out + '_1_domain_notrepeated', 'w')
    Domain_file_modif_complex = open(Domain_file_out + '_1_domain_complex_modif', 'w')
    Domain_file_repeated = open(Domain_file_out + '_1_domain_repeated', 'w')
    #
    
    #store the pair in the appropriate final group
    PairOrtho = open(PairOrtho_1domain_in, 'r')
    for pair_ortho in PairOrtho:
        domain_sp1 = spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]]
        domain_sp2 = spec_domain[pair_ortho.split('\t')[2].rstrip('0123456789')][pair_ortho.split('\t')[2]]
        
        if set(domain_sp1) == set(domain_sp2):
            Domain_file_repeated.write(pair_ortho)
        elif len(set(domain_sp1)) + 1 == len(set(domain_sp2)):
            not_all_present = False
            for domain_considered in set(domain_sp1):
                if domain_considered not in set(domain_sp2):
                    not_all_present = True
                    break
            if not_all_present:
                Domain_file_modif_complex.write(pair_ortho)
            else:
                Domain_file_notrepeated.write(pair_ortho)
        elif len(set(domain_sp2)) + 1 == len(set(domain_sp1)):
            not_all_present = False
            for domain_considered in set(domain_sp2):
                if domain_considered not in set(domain_sp1):
                    not_all_present = True
                    break
            if not_all_present:
                Domain_file_modif_complex.write(pair_ortho)
            else:
                Domain_file_notrepeated.write(pair_ortho)
        else:
            Domain_file_modif_complex.write(pair_ortho)
                
    PairOrtho.close()
    Domain_file_notrepeated.close()
    Domain_file_modif_complex.close()
    Domain_file_repeated.close()
    #
    
#Run the function / WINDOWS
#In different mammals
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO'], 
                                    central_mammals = 'HUMAN'):
    for considered_mammals in list_mammals:
        extract_ortholog_repeated_domain(PairOrtho_1domain_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s_1_domain_modif' % (central_mammals, considered_mammals),
                                         Domain_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s_domain' % (central_mammals, considered_mammals),
                                         Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s' % (central_mammals, considered_mammals),
                                         central_species = central_mammals,
                                         considered_species = considered_mammals)    
#
run_function_in_mammals_dataset()