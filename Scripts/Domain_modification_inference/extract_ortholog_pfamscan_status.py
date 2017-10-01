# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 15:25:57 2017

@author: Claivaz
extraction of pair which have at least one domain inferred by pfamscan (both/partial/none), usefull for analysis of domain architecture inference
"""

def extract_ortholog_pfamscan_status(PairOrtho_in, 
                                     Domain_file_in, 
                                     Domain_file_out,
                                     central_species,
                                     considered_species,
                                     parameter = '1:1'):
    
    PairOrtho = open(PairOrtho_in, 'r')
    
    #create gene dictionnary
    spec_domain = {}
    spec_domain[central_species] = {}
    spec_domain[considered_species] = {}
    
    for pair_ortho in PairOrtho:
        if parameter in pair_ortho:
            spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]] = []
            spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]] = []
    
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
    Domain_file_partial = open(Domain_file_out + '_partial', 'w')
    Domain_file_none = open(Domain_file_out + '_none', 'w')
    Domain_file_both = open(Domain_file_out + '_both', 'w')
    #
    
    #check if the considered pair have at least one domain inferred by pfamscan and store the pair in the appropriate final group
    PairOrtho = open(PairOrtho_in, 'r')
    for pair_ortho in PairOrtho:    
        if parameter in pair_ortho:
            if central_species in pair_ortho.split('\t')[0] and considered_species in pair_ortho.split('\t')[1] or central_species in pair_ortho.split('\t')[1] and considered_species in pair_ortho.split('\t')[0]:
                if len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]]) > 0 and len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]]) > 0:                
                    Domain_file_both.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])) +
                                              '\t' + pair_ortho.split('\t')[1] + '\tboth\n')
                elif len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]]) > 0 or len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]]) > 0:
                    Domain_file_partial.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])) +
                                              '\t' + pair_ortho.split('\t')[1] + '\tpartial\n')
                else:
                    Domain_file_none.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])) +
                                              '\t' + pair_ortho.split('\t')[1] + '\tnone\n')
                    
    PairOrtho.close()
    Domain_file_none.close()
    Domain_file_both.close()
    Domain_file_partial.close()
    #
    
#Run the function / WINDOWS
#In different mammals
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX'], 
                                    central_mammals = 'HUMAN'):
    for considered_mammals in list_mammals:
        extract_ortholog_pfamscan_status(PairOrtho_in = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_%s_%s.txt' % (central_mammals, considered_mammals),
                                         Domain_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s_domain' % (central_mammals, considered_mammals),
                                         Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s' % (central_mammals, considered_mammals),
                                         central_species = central_mammals,
                                         considered_species = considered_mammals)    
#
run_function_in_mammals_dataset()
