# -*- coding: utf-8 -*-
"""
171009

@author: Claivaz
analysis of domain modification inference amongst mammals paralog
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
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN']):
    dict_results = {}
    path_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/'

    for considered_mammals in list_mammals:
        
        dict_results[considered_mammals] = {}
        
        #Unique ortholog call from OMA: OMA
        dict_results[considered_mammals]['OMA'] = number_pair_in_file('D:/UNIL/Master/Master_Project/Data/OMA/pairwise_paralog_%s' % (considered_mammals) )
        
        #Domain inferred by pfamscan (none/both/partial): pfamscan
        dict_results[considered_mammals]['pfamscan'] = {}
        
        dict_results[considered_mammals]['pfamscan']['partial'] = number_pair_in_file(path_file + 'paralog_%s_partial' % (considered_mammals))
        dict_results[considered_mammals]['pfamscan']['both'] = number_pair_in_file(path_file + 'paralog_%s_both' % (considered_mammals))
        dict_results[considered_mammals]['pfamscan']['none'] = number_pair_in_file(path_file + 'paralog_%s_none' % (considered_mammals))
            
        #Domain modification group (no modif/complex_modif/1_domain_modif): group
        dict_results[considered_mammals]['group'] = {}
        
        dict_results[considered_mammals]['group']['nomodif'] = number_pair_in_file(path_file + 'paralog_%s_domain_nomodif' % (considered_mammals))
        dict_results[considered_mammals]['group']['complex_modif'] = number_pair_in_file(path_file + 'paralog_%s_domain_complex_modif' % (considered_mammals))
            
        #Domain repetition (1_domain_notrepeated/1_domain_repeated/1_domain_complex_modif): repeat
    
        dict_results[considered_mammals]['repeat'] = {}
        dict_results[considered_mammals]['repeat']['1_domain_notrepeated'] = number_pair_in_file(path_file + 'paralog_%s_1_domain_notrepeated' % (considered_mammals))
        dict_results[considered_mammals]['repeat']['1_domain_repeated'] = number_pair_in_file(path_file + 'paralog_%s_1_domain_repeated' % (considered_mammals))
            
        #Intersect domainDIFF and 1_domain_notrepeated: final
        dict_results[considered_mammals]['final'] = {} 
        
        dict_results[considered_mammals]['final']['f-1'] = number_pair_in_file(path_file + 'putative_paralog_%s_domain_loss' % (considered_mammals), reg_exp = 'f-1')
        dict_results[considered_mammals]['final']['b-1'] = number_pair_in_file(path_file + 'putative_paralog_%s_domain_loss' % (considered_mammals), reg_exp = 'b-1')
        dict_results[considered_mammals]['final']['int-1'] = number_pair_in_file(path_file + 'putative_paralog_%s_domain_loss' % (considered_mammals), reg_exp = 'int-1')
    
    return(dict_results)
#
run_function_in_mammals_dataset()
