# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 11:40:17 2017

@author: Claivaz

Inference of loss domain on specific tissue expression, domain family (clan) enrichment
PfamA clan file can be downloaded from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
"""

def clan_inference(PairOrtho_in_path, 
                   Domain_file_in_path,
                   second_species,
                   considered_species,
                   clan_file,
                   TspecOrtho_path,
                   parameter = '1:1'):

    #create gene dictionnary
    try:
        PairOrtho = open('%s_%s_%s' % (PairOrtho_in_path, considered_species, second_species), 'r')
    
    except:
        PairOrtho = open('%s_%s_%s' % (PairOrtho_in_path, second_species, considered_species), 'r')
    
    spec_domain = {}
    spec_domain[second_species] = {}
    spec_domain[considered_species] = {}
    
    for pair_ortho in PairOrtho:
        
        if parameter in pair_ortho:
            spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]] = []
            spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]] = []
    
    PairOrtho.close()
    #

    #store domain into dictionnary
    for specie_tmp in [considered_species, second_species]:   
        Domain_file = open('%s%s_domain' % (Domain_file_in_path, specie_tmp), 'r')
        
        for ortho_domain_line in Domain_file: 
            if ortho_domain_line.replace(' ','\t').split('\t')[0] in spec_domain[second_species].keys() or ortho_domain_line.replace(' ','\t').split('\t')[0] in spec_domain[considered_species].keys():
                ortho_domain_line = ortho_domain_line.replace(' ', '\t')
                
                while '\t\t' in ortho_domain_line:
                    ortho_domain_line = ortho_domain_line.replace('\t\t', '\t')
                
                spec_domain[ortho_domain_line.split('\t')[0].rstrip('0123456789')][ortho_domain_line.split('\t')[0]].append(ortho_domain_line.split('\t')[6])
        
        Domain_file.close()
    #
    
    #store clan into dictionnary
    clan_dict = {}
    clan_open = open(clan_file, 'r')
    
    for clan_line in clan_open:
        clan_dict[clan_line.split('\t')[3]] = clan_line.split('\t')[1]
        
    clan_open.close()
    #
    
    #open output of Tspec_inference.R and store all the tissue available
    try:
        TspecOrtho = open('%s%s_%s_ortho_notfemale_dataset' % (TspecOrtho_path, considered_species, second_species), 'r')
        
    except:
        TspecOrtho = open('%s%s_%s_ortho_notfemale_dataset' % (TspecOrtho_path, second_species, considered_species), 'r')
    
    TspecOrtho.readline()
    
    tissue_list = []
    for Tspec_line in TspecOrtho:
        if Tspec_line.replace('"', '').split(',')[8] not in tissue_list:
            tissue_list.append(Tspec_line.replace('"', '').split(',')[8])
            
        if Tspec_line.replace('"', '').split(',')[11] not in tissue_list:
            tissue_list.append(Tspec_line.replace('"', '').split(',')[11])
    
    #create dictionnary with tissue
    tissue_dict = {}
    
    for tissue_tmp in tissue_list:
        tissue_dict[tissue_tmp] = []
        
    #open output of Tspec_inference.R    
    try:
        TspecOrtho = open('%s%s_%s_ortho_notfemale_dataset' % (TspecOrtho_path, considered_species, second_species), 'r')
        
    except:
        TspecOrtho = open('%s%s_%s_ortho_notfemale_dataset' % (TspecOrtho_path, second_species, considered_species), 'r')
    
    TspecOrtho.readline()
    
    for Tspec_line in TspecOrtho:
        #consider only pairs with a shift from ubiquitous to specific which can be due to domain loss/gain
        if 'control' not in Tspec_line and 'ubiquitous' in Tspec_line and 'specific' in Tspec_line:
            if clan_dict[list(set(spec_domain[Tspec_line.replace('"', '').split(',')[1].rstrip('0123456789')][Tspec_line.replace('"', '').split(',')[1]]).symmetric_difference(set(spec_domain[Tspec_line.replace('"', '').split(',')[3].rstrip('0123456789')][Tspec_line.replace('"', '').split(',')[3]])))[0]] != '':
                if Tspec_line.replace('"', '').split(',')[11] == 'ubiquitous':
                    tissue_dict[Tspec_line.replace('"', '').split(',')[8]].append(clan_dict[list(set(spec_domain[Tspec_line.replace('"', '').split(',')[1].rstrip('0123456789')][Tspec_line.replace('"', '').split(',')[1]]).symmetric_difference(set(spec_domain[Tspec_line.replace('"', '').split(',')[3].rstrip('0123456789')][Tspec_line.replace('"', '').split(',')[3]])))[0]])
                
                else:
                    tissue_dict[Tspec_line.replace('"', '').split(',')[11]].append(clan_dict[list(set(spec_domain[Tspec_line.replace('"', '').split(',')[1].rstrip('0123456789')][Tspec_line.replace('"', '').split(',')[1]]).symmetric_difference(set(spec_domain[Tspec_line.replace('"', '').split(',')[3].rstrip('0123456789')][Tspec_line.replace('"', '').split(',')[3]])))[0]])
        
    TspecOrtho.close()
    #
    
    from collections import Counter
    
    for key_tmp in tissue_dict.keys():
        if key_tmp != 'ubiquitous':
            freqs_tmp = Counter(tissue_dict[key_tmp])
            print('%s:\n%s' % (key_tmp, freqs_tmp))
    
    freqs_tmp = []                
    for key_tmp in tissue_dict.keys():
        if key_tmp != 'ubiquitous':
            freqs_tmp.append(tissue_dict[key_tmp])    
    from itertools import chain
    print(Counter(list(chain.from_iterable(freqs_tmp))))
    
    
#Run in Windows for the whole datasets   
mammals_done = []
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN']):
    for considered_mammals in list_mammals:
        mammals_done.append(considered_mammals)
        for second_mammals in list_mammals:
            if second_mammals not in mammals_done:
                clan_inference(PairOrtho_in_path = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog',
                               Domain_file_in_path = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/',
                               TspecOrtho_path = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/',
                               second_species = second_mammals,
                               considered_species = considered_mammals,
                               clan_file = 'D:/UNIL/Master/Master_Project/Data/pfam/Pfam-A.clans.tsv')
#Run in windows just one pairwise
clan_inference(PairOrtho_in_path = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog',
               Domain_file_in_path = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/',
               TspecOrtho_path = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/',
               second_species = 'MOUSE',
               considered_species = 'HUMAN',
               clan_file = 'D:/UNIL/Master/Master_Project/Data/pfam/Pfam-A.clans.tsv')
