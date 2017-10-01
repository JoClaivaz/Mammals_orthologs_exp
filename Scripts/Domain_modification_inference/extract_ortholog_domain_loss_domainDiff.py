# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170522
Extraction of the unique ortholog pair, considering that the loss it's not part of repetition domains and inferred by domainDIFF
"""
def extract_ortholog_domain_loss_domainDiff (domain_loss_in, domain_modification_in, output_file):
    
    domain_loss = open(domain_loss_in, 'r')
    pair_out = open(output_file, 'w')
    domain_modification = open(domain_modification_in, 'r')
    
    dict_domainloss = {}    
    for domain_loss_line in domain_loss:
        dict_domainloss[domain_loss_line.split('\t')[0]] = domain_loss_line.split('\t')[2]
        dict_domainloss[domain_loss_line.split('\t')[2]] = domain_loss_line.split('\t')[0]
    
    domain_loss.close()
    
    for domain_modification_line in domain_modification:
        try:
            dict_domainloss[domain_modification_line.split(' ')[0]]
            
            if dict_domainloss[domain_modification_line.split(' ')[0]] == domain_modification_line.split(' ')[2]:
                pair_out.write(domain_modification_line.replace(' ', '\t'))
        except:
            pass
        
    domain_modification.close()
    pair_out.close()
        
#Run the function / WINDOWS
#In different mammals
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX'], 
                                    central_mammals = 'HUMAN'):
    for considered_mammals in list_mammals:
        extract_ortholog_domain_loss_domainDiff(domain_loss_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s_1_domain_notrepeated' % (central_mammals, considered_mammals), 
                                                domain_modification_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s_domain_modifications' % (central_mammals, considered_mammals), 
                                                output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/final_pair_%s_%s_domain_loss' % (central_mammals, considered_mammals))    
#
run_function_in_mammals_dataset()
