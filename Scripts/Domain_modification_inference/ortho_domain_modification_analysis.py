# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 15:51:25 2017

@author: Claivaz
analysis of domain modification inference amongst mammals ortholog
"""

def number_pair_in_file(file_in, reg_exp = False):
    file_open = open(file_in, 'r')
    count_line = 0
    for file_line in file_open:
        if not reg_exp: 
            count_line += 1
        elif reg_exp and reg_exp in file_line:
            count_line += 1
            
    file_open.close()
    
    return(count_line)


#Run the function / WINDOWS
#In different mammals
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO'], 
                                    central_mammals = 'HUMAN'):
    dict_results = {}
    path_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/'

    for considered_mammals in list_mammals:
        
        dict_results[considered_mammals] = {}
        
        #Unique ortholog call from OMA: OMA
        try:
            dict_results[considered_mammals]['OMA'] = number_pair_in_file('D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_%s_%s' % (central_mammals, considered_mammals), reg_exp = '1:1' )
        
        except:
            dict_results[considered_mammals]['OMA'] = number_pair_in_file('D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_%s_%s' % (considered_mammals, central_mammals), reg_exp = '1:1' )
            
        #Domain inferred by pfamscan (none/both/partial): pfamscan
        dict_results[considered_mammals]['pfamscan'] = {}
        
        try: 
            dict_results[considered_mammals]['pfamscan']['partial'] = number_pair_in_file(path_file + 'ortholog_%s_%s_partial' % (central_mammals, considered_mammals))
            dict_results[considered_mammals]['pfamscan']['both'] = number_pair_in_file(path_file + 'ortholog_%s_%s_both' % (central_mammals, considered_mammals))
            dict_results[considered_mammals]['pfamscan']['none'] = number_pair_in_file(path_file + 'ortholog_%s_%s_none' % (central_mammals, considered_mammals))
            
        except:
            dict_results[considered_mammals]['pfamscan']['partial'] = number_pair_in_file(path_file + 'ortholog_%s_%s_partial' % (considered_mammals, central_mammals))
            dict_results[considered_mammals]['pfamscan']['both'] = number_pair_in_file(path_file + 'ortholog_%s_%s_both' % (considered_mammals, central_mammals))
            dict_results[considered_mammals]['pfamscan']['none'] = number_pair_in_file(path_file + 'ortholog_%s_%s_none' % (considered_mammals, central_mammals))
            
        #Domain modification group (no modif/complex_modif/1_domain_modif): group
        dict_results[considered_mammals]['group'] = {}
        
        try:
            dict_results[considered_mammals]['group']['nomodif'] = number_pair_in_file(path_file + 'ortholog_%s_%s_domain_nomodif' % (central_mammals, considered_mammals))
            dict_results[considered_mammals]['group']['complex_modif'] = number_pair_in_file(path_file + 'ortholog_%s_%s_domain_complex_modif' % (central_mammals, considered_mammals))
            
        except:
            dict_results[considered_mammals]['group']['nomodif'] = number_pair_in_file(path_file + 'ortholog_%s_%s_domain_nomodif' % (considered_mammals, central_mammals))
            dict_results[considered_mammals]['group']['complex_modif'] = number_pair_in_file(path_file + 'ortholog_%s_%s_domain_complex_modif' % (considered_mammals, central_mammals))
            
        #Domain repetition (1_domain_notrepeated/1_domain_repeated/1_domain_complex_modif): repeat
        try:
            dict_results[considered_mammals]['repeat'] = {}
            dict_results[considered_mammals]['repeat']['1_domain_notrepeated'] = number_pair_in_file(path_file + 'ortholog_%s_%s_1_domain_notrepeated' % (central_mammals, considered_mammals))
            dict_results[considered_mammals]['repeat']['1_domain_repeated'] = number_pair_in_file(path_file + 'ortholog_%s_%s_1_domain_repeated' % (central_mammals, considered_mammals))
        
        except:
            dict_results[considered_mammals]['repeat'] = {}
            dict_results[considered_mammals]['repeat']['1_domain_notrepeated'] = number_pair_in_file(path_file + 'ortholog_%s_%s_1_domain_notrepeated' % (considered_mammals, central_mammals))
            dict_results[considered_mammals]['repeat']['1_domain_repeated'] = number_pair_in_file(path_file + 'ortholog_%s_%s_1_domain_repeated' % (considered_mammals, central_mammals))
            
        #Intersect domainDIFF and 1_domain_notrepeated: final
        dict_results[considered_mammals]['final'] = {} 
        
        try:
            dict_results[considered_mammals]['final']['f-1'] = number_pair_in_file(path_file + 'putative_ortholog_%s_%s_domain_loss' % (central_mammals, considered_mammals), reg_exp = 'f-1')
            dict_results[considered_mammals]['final']['b-1'] = number_pair_in_file(path_file + 'putative_ortholog_%s_%s_domain_loss' % (central_mammals, considered_mammals), reg_exp = 'b-1')
            dict_results[considered_mammals]['final']['int-1'] = number_pair_in_file(path_file + 'putative_ortholog_%s_%s_domain_loss' % (central_mammals, considered_mammals), reg_exp = 'int-1')
        
        except:
            dict_results[considered_mammals]['final']['f-1'] = number_pair_in_file(path_file + 'putative_ortholog_%s_%s_domain_loss' % (considered_mammals, central_mammals), reg_exp = 'f-1')
            dict_results[considered_mammals]['final']['b-1'] = number_pair_in_file(path_file + 'putative_ortholog_%s_%s_domain_loss' % (considered_mammals, central_mammals), reg_exp = 'b-1')
            dict_results[considered_mammals]['final']['int-1'] = number_pair_in_file(path_file + 'putative_ortholog_%s_%s_domain_loss' % (considered_mammals, central_mammals), reg_exp = 'int-1')
        
    return(dict_results)
#
run_function_in_mammals_dataset()
