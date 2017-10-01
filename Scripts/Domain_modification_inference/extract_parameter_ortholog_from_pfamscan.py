# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170520
Extraction of pfamscan data, in function of pairwise_OMA and add pair_number (usefull for next analysis steps)
"""
def extract_parameter_ortholog_from_pfamscan(pairwise_oma_in, pfamscan_in, all_ortholog_out, parameter = '1:1'):
    
    pairwise_oma = open(pairwise_oma_in, 'r')
    domain_out = open(all_ortholog_out, 'w')
    
    pair_number = 0
    for pair in pairwise_oma:
        if parameter in pair.split('\t')[2] :
            pair_number += 1
            pfamscan_out = open(pfamscan_in, 'r')
            for domain_line in pfamscan_out:
                
                if pair.split('\t')[0] in domain_line:
                    domain_line = domain_line.replace(' ', '\t')
                    while '\t\t' in domain_line:
                        domain_line = domain_line.replace('\t\t', '\t')
                    domain_out.write(str(pair_number) + '\t' + domain_line.split('\t')[0] + '\t' +
                                     domain_line.split('\t')[5] + '\t' + domain_line.split('\t')[8] + '\t' +
                                     domain_line.split('\t')[9] + '\t' + domain_line.split('\t')[10] + '\n')
                
                elif pair.split('\t')[1] in domain_line:
                    domain_line = domain_line.replace(' ', '\t')
                    while '\t\t' in domain_line:
                        domain_line = domain_line.replace('\t\t', '\t')
                    domain_out.write(str(pair_number) + '\t' + domain_line.split('\t')[0] + '\t' +
                                     domain_line.split('\t')[5] + '\t' + domain_line.split('\t')[8] + '\t' +
                                     domain_line.split('\t')[9] + '\t' + domain_line.split('\t')[10] + '\n')
            pfamscan_out.close()
    pairwise_oma.close()
    domain_out.close()            
    
#Run the function / WINDOWS
#In different mammals
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX'], 
                                    central_mammals = 'HUMAN'):
    for considered_mammals in list_mammals:
        extract_parameter_ortholog_from_pfamscan(pairwise_oma_in = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_%s_%s.txt' % (central_mammals, considered_mammals), 
                                                 pfamscan_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_%s_%s_domain' % (central_mammals, considered_mammals),
                                                 all_ortholog_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/unique_ortholog_%s_%s_domain_parsed' % (central_mammals, considered_mammals))    
#
run_function_in_mammals_dataset()