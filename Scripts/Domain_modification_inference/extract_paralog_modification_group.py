# -*- coding: utf-8 -*-
"""
171001, last modification 171009

@author: Claivaz

extraction of pair which have at least one domain inferred by pfamscan (both/partial/none)
extraction of pair which have only one domain modification and no modification (control group)
extraction of pair with domain modification not involved in domain repetition
Inference of the position of the loss (f-1, b-1, int-1; N-terminal, C-terminal, internal respectively)

Command line
python extract_paralog_modification_group.py list_mammals PairPara_in_path_with_regexp Domain_file_in_path Domain_file_out_with_regexp
"""

def extract_paralog_modification_group(PairPara_in_path, 
                                       Domain_file_in_path, 
                                       Domain_file_out,
                                       considered_species):
    
    
    PairPara = open('%s_%s' % (PairPara_in_path, considered_species), 'r')
    
    #create gene dictionnary
    spec_domain = {}
    
    for pair_para in PairPara:
        spec_domain[pair_para.split('\t')[0]] = []
        spec_domain[pair_para.split('\t')[1]] = []
    
    PairPara.close()
    #
    
    #store domain into dictionnary
    Domain_file = open('%s%s_domain' % (Domain_file_in_path, considered_species), 'r')
    
    for para_domain_line in Domain_file:
        
        if para_domain_line.replace(' ','\t').split('\t')[0] in spec_domain.keys():
            para_domain_line = para_domain_line.replace(' ', '\t')
            
            while '\t\t' in para_domain_line:
                para_domain_line = para_domain_line.replace('\t\t', '\t')
            
            spec_domain[para_domain_line.split('\t')[0]].append(para_domain_line.split('\t')[6])
    
    Domain_file.close()
    #   
        
    #open file for different group modification    
    Domain_file_partial = open(Domain_file_out + '_partial', 'w')
    Domain_file_none = open(Domain_file_out + '_none', 'w')
    Domain_file_both = open(Domain_file_out + '_both', 'w')
    Domain_file_nomodif = open(Domain_file_out + '_domain_nomodif', 'w')
    Domain_file_modif_complex = open(Domain_file_out + '_domain_complex_modif', 'w')
    Domain_file_modif_1 = open(Domain_file_out + '_1_domain_modif', 'w')
    Domain_file_notrepeated = open(Domain_file_out + '_1_domain_notrepeated', 'w')
    Domain_file_repeated = open(Domain_file_out + '_1_domain_repeated', 'w')
    Domain_file_final = open('%sputative_paralog_%s_domain_loss' % (Domain_file_out.split('paralog_')[0], considered_species), 'w')
    #
    
    #sorting paralog pair
    PairPara = open('%s_%s' % (PairPara_in_path, considered_species), 'r')
    
    for pair_para in PairPara:    
        
        #check if the considered pair have at least one domain inferred by pfamscan and store the pair in the appropriate final group
        if len(spec_domain[pair_para.split('\t')[0]]) > 0 and len(spec_domain[pair_para.split('\t')[1]]) > 0:                
            Domain_file_both.write(pair_para.split('\t')[0] + '\t' + str(len(spec_domain[pair_para.split('\t')[1]])) +
                                      '\t' + pair_para.split('\t')[1] + '\tboth\n')
            
            #check if the considered pair have unmodified domain architecture or one domain lost and store the pair in appropriate file
            if spec_domain[pair_para.split('\t')[0]] == spec_domain[pair_para.split('\t')[1]]:
                Domain_file_nomodif.write(pair_para.split('\t')[0] + '\t' +  pair_para.split('\t')[3] + '\t' +
                                          pair_para.split('\t')[1] + '\t' + pair_para.split('\t')[4].replace('\n', '') + '\tnomodif\n')
            elif len(spec_domain[pair_para.split('\t')[0]]) == len(spec_domain[pair_para.split('\t')[1]]):
                Domain_file_modif_complex.write(pair_para.split('\t')[0] + '\t' +  pair_para.split('\t')[3] + '\t' +
                                                pair_para.split('\t')[1] + '\t' + pair_para.split('\t')[4].replace('\n', '') + '\tcomplex_modif\n')
            else:
                if len(spec_domain[pair_para.split('\t')[0]]) + 1 == len(spec_domain[pair_para.split('\t')[1]]) or len(spec_domain[pair_para.split('\t')[0]]) - 1 == len(spec_domain[pair_para.split('\t')[1]]):
                    Domain_file_modif_1.write(pair_para.split('\t')[0] + '\t' +  pair_para.split('\t')[3] + '\t' +
                                              pair_para.split('\t')[1] + '\t' + pair_para.split('\t')[4].replace('\n', '') + '\tmodif\n')

                    #check if the loss domain is repeated into the protein or not                    
                    domain_sp1 = spec_domain[pair_para.split('\t')[0]]
                    domain_sp2 = spec_domain[pair_para.split('\t')[1]]                    
                    
                    if set(domain_sp1) == set(domain_sp2):
                        Domain_file_repeated.write(pair_para)
                    
                    elif len(set(domain_sp1)) + 1 == len(set(domain_sp2)):
                        not_all_present = False
                        
                        for domain_considered in set(domain_sp1):
                            
                            if domain_considered not in set(domain_sp2):
                                not_all_present = True
                                break
                        
                        if not_all_present:
                            Domain_file_modif_complex.write(pair_para)
                        
                        else:
                            Domain_file_notrepeated.write(pair_para)
                            
                            #Inference of the position of the loss
                            if domain_sp1[1:] == domain_sp2 or domain_sp2[1:] == domain_sp1:
                                Domain_file_final.write('%s\tf-1\t%s\t%s\t%s\n' % (pair_para.split('\t')[0], pair_para.split('\t')[1], pair_para.split('\t')[3], pair_para.split('\t')[4].replace('\n', '')))
                                
                            elif domain_sp1[:-1] == domain_sp2 or domain_sp2[:-1] == domain_sp1:
                                Domain_file_final.write('%s\tb-1\t%s\t%s\t%s\n' % (pair_para.split('\t')[0], pair_para.split('\t')[1], pair_para.split('\t')[3], pair_para.split('\t')[4].replace('\n', '')))
                            
                            else:
                                Domain_file_final.write('%s\tint-1\t%s\t%s\t%s\n' % (pair_para.split('\t')[0], pair_para.split('\t')[1], pair_para.split('\t')[3], pair_para.split('\t')[4].replace('\n', '')))
                                                        
                    elif len(set(domain_sp2)) + 1 == len(set(domain_sp1)):
                        not_all_present = False
                        
                        for domain_considered in set(domain_sp2):
                            
                            if domain_considered not in set(domain_sp1):
                                not_all_present = True
                                break
                        
                        if not_all_present:
                            Domain_file_modif_complex.write(pair_para)
                        
                        else:
                            Domain_file_notrepeated.write(pair_para)
                            
                            #Inference of the position of the loss
                            if domain_sp1[1:] == domain_sp2 or domain_sp2[1:] == domain_sp1:
                                Domain_file_final.write('%s\tf-1\t%s\t%s\t%s\n' % (pair_para.split('\t')[0], pair_para.split('\t')[1], pair_para.split('\t')[3], pair_para.split('\t')[4].replace('\n', '')))
                                
                            elif domain_sp1[:-1] == domain_sp2 or domain_sp2[:-1] == domain_sp1:
                                Domain_file_final.write('%s\tb-1\t%s\t%s\t%s\n' % (pair_para.split('\t')[0], pair_para.split('\t')[1], pair_para.split('\t')[3], pair_para.split('\t')[4].replace('\n', '')))
                            
                            else:
                                Domain_file_final.write('%s\tint-1\t%s\t%s\t%s\n' % (pair_para.split('\t')[0], pair_para.split('\t')[1], pair_para.split('\t')[3], pair_para.split('\t')[4].replace('\n', '')))
                            
                    else:
                        Domain_file_modif_complex.write(pair_para)
                                        
                    
                else:
                    Domain_file_modif_complex.write(pair_para.split('\t')[0] + '\t' + str(len(spec_domain[pair_para.split('\t')[0]])) +
                                                    '\t' + pair_para.split('\t')[1] + '\tcomplex_modif\n')
            
            
            
        elif len(spec_domain[pair_para.split('\t')[0]]) > 0 or len(spec_domain[pair_para.split('\t')[1]]) > 0:
            Domain_file_partial.write(pair_para.split('\t')[0] + '\t' + str(len(spec_domain[pair_para.split('\t')[1]])) +
                                      '\t' + pair_para.split('\t')[1] + '\tpartial\n')
        else:
            Domain_file_none.write(pair_para.split('\t')[0] + '\t' + str(len(spec_domain[pair_para.split('\t')[1]])) +
                                      '\t' + pair_para.split('\t')[1] + '\tnone\n')
            
    PairPara.close()
    Domain_file_none.close()
    Domain_file_both.close()
    Domain_file_partial.close()
    Domain_file_nomodif.close()
    Domain_file_modif_complex.close()
    Domain_file_modif_1.close()
    Domain_file_notrepeated.close()
    Domain_file_repeated.close()
    Domain_file_final.close()
    #
    
#Run the function / from command line
#python extract_paralog_modification_group.py list_mammals PairPara_in_path_with_regexp Domain_file_in_path Domain_file_out_with_regexp
import sys
def run_function_in_mammals_dataset(list_mammals = sys.argv[1]):
    
    for considered_mammals in list_mammals:
        extract_paralog_modification_group(PairPara_in_path = sys.argv[2],
                                           Domain_file_in_path = sys.argv[3],
                                           Domain_file_out = '%s%s' % (sys.argv[4], considered_mammals),
                                           considered_species = considered_mammals)
''' Run from python
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN']):
    for considered_mammals in list_mammals:
        extract_paralog_modification_group(PairPara_in_path = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_paralog',
                                           Domain_file_in_path = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/',
                                           Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_%s' % (considered_mammals),
                                           considered_species = considered_mammals)
'''
#
run_function_in_mammals_dataset()
