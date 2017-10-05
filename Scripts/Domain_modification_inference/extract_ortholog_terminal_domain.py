# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:49:25 2017

@author: Claivaz
Extract pair with domain modification in the terminal part of the protein
"""


def extract_ortholog_terminal_domain(PairOrtho_1domain_notrepeated, 
                                     Domain_file_in, 
                                     Domain_file_out,
                                     central_species,
                                     considered_species):
    
    PairOrtho = open(PairOrtho_1domain_notrepeated, 'r')
    
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
    Domain_file_terminal = open(Domain_file_out, 'w')
    #
    
    #store the pair in the appropriate final group
    PairOrtho = open(PairOrtho_1domain_notrepeated, 'r')
    for pair_ortho in PairOrtho:
        domain_sp1 = spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]]
        domain_sp2 = spec_domain[pair_ortho.split('\t')[2].rstrip('0123456789')][pair_ortho.split('\t')[2]]
        
        if domain_sp1[1:] == domain_sp2 or domain_sp2[1:] == domain_sp1:
            Domain_file_terminal.write('%s\t%s\t%s\tf-1\n' % (pair_ortho.split('\t')[0], pair_ortho.split('\t')[1], pair_ortho.split('\t')[2]))
            
        elif domain_sp1[:-1] == domain_sp2 or domain_sp2[:-1] == domain_sp1:
            Domain_file_terminal.write('%s\t%s\t%s\tb-1\n' % (pair_ortho.split('\t')[0], pair_ortho.split('\t')[1], pair_ortho.split('\t')[2]))
            
                
    PairOrtho.close()
    Domain_file_terminal.close()
    #

#Run the function / WINDOWS
#In different mammals
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO'], 
                                    central_mammals = 'HUMAN'):
    for considered_mammals in list_mammals:
        extract_ortholog_terminal_domain(PairOrtho_1domain_notrepeated = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s_1_domain_notrepeated' % (central_mammals, considered_mammals),
                                         Domain_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s_domain' % (central_mammals, considered_mammals),
                                         Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/final_pair_%s_%s_domain_loss' % (central_mammals, considered_mammals),
                                         central_species = central_mammals,
                                         considered_species = considered_mammals)    
#
run_function_in_mammals_dataset()