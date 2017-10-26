# -*- coding: utf-8 -*-
"""
171026

@author: Claivaz

extraction of pair which have at least one domain inferred by pfamscan (both/partial/none)
extraction of pair which have only one domain modification and no modification (control group)
extraction of pair with domain modification not involved in domain repetition
Inference of the position of the loss (f-1, b-1, int-1; N-terminal, C-terminal, internal respectively)
The parameter '1:1' indicate unique ortholog, if the consideration is on whole ortholog parameter should be set as ':'
"""

def extract_paraortho_modification_group(PairOrtho_in_path, 
                                         Domain_file_in_path, 
                                         Domain_file_out,
                                         second_species,
                                         considered_species,
                                         with_testis,
                                         parameter = ':'):
    
    if with_testis:
        try:
            PairOrtho = open('%s_%s_%s_withtestis' % (PairOrtho_in_path, considered_species, second_species), 'r')
        
        except:
            PairOrtho = open('%s_%s_%s_withtestis' % (PairOrtho_in_path, second_species, considered_species), 'r')
    else:
        try:
            PairOrtho = open('%s_%s_%s_withouttestis' % (PairOrtho_in_path, considered_species, second_species), 'r')
        
        except:
            PairOrtho = open('%s_%s_%s_withouttestis' % (PairOrtho_in_path, second_species, considered_species), 'r')
        
    #create gene dictionnary
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
        
    #open file for different group modification    
    if with_testis:
        Domain_file_partial = open(Domain_file_out + '_partial_withtestis', 'w')
        Domain_file_none = open(Domain_file_out + '_none_withtestis', 'w')
        Domain_file_both = open(Domain_file_out + '_both_withtestis', 'w')
        Domain_file_nomodif = open(Domain_file_out + '_domain_nomodif_withtestis', 'w')
        Domain_file_modif_complex = open(Domain_file_out + '_domain_complex_modif_withtestis', 'w')
        Domain_file_modif_1 = open(Domain_file_out + '_1_domain_modif_withtestis', 'w')
        Domain_file_notrepeated = open(Domain_file_out + '_1_domain_notrepeated_withtestis', 'w')
        Domain_file_repeated = open(Domain_file_out + '_1_domain_repeated_withtestis', 'w')
        Domain_file_final = open('%sputative_paraortho_%s_%s_domain_loss_withtestis' % (Domain_file_out.split('paraortho_')[0], second_species, considered_species), 'w')
    else:
        Domain_file_partial = open(Domain_file_out + '_partial_withouttestis', 'w')
        Domain_file_none = open(Domain_file_out + '_none_withouttestis', 'w')
        Domain_file_both = open(Domain_file_out + '_both_withouttestis', 'w')
        Domain_file_nomodif = open(Domain_file_out + '_domain_nomodif_withouttestis', 'w')
        Domain_file_modif_complex = open(Domain_file_out + '_domain_complex_modif_withouttestis', 'w')
        Domain_file_modif_1 = open(Domain_file_out + '_1_domain_modif_withouttestis', 'w')
        Domain_file_notrepeated = open(Domain_file_out + '_1_domain_notrepeated_withouttestis', 'w')
        Domain_file_repeated = open(Domain_file_out + '_1_domain_repeated_withouttestis', 'w')
        Domain_file_final = open('%sputative_paraortho_%s_%s_domain_loss_withouttestis' % (Domain_file_out.split('paraortho_')[0], second_species, considered_species), 'w')
    #
    
    #sorting ortholog pair
    if with_testis:
        try:
            PairOrtho = open('%s_%s_%s_withtestis' % (PairOrtho_in_path, considered_species, second_species), 'r')
        
        except:
            PairOrtho = open('%s_%s_%s_withtestis' % (PairOrtho_in_path, second_species, considered_species), 'r')
    else:
        try:
            PairOrtho = open('%s_%s_%s_withouttestis' % (PairOrtho_in_path, considered_species, second_species), 'r')
        
        except:
            PairOrtho = open('%s_%s_%s_withouttestis' % (PairOrtho_in_path, second_species, considered_species), 'r')
    
    for pair_ortho in PairOrtho:    
        
        if parameter in pair_ortho:
            
            if second_species in pair_ortho.split('\t')[0] and considered_species in pair_ortho.split('\t')[1] or second_species in pair_ortho.split('\t')[1] and considered_species in pair_ortho.split('\t')[0]:
                
                #check if the considered pair have at least one domain inferred by pfamscan and store the pair in the appropriate final group
                if len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]]) > 0 and len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]]) > 0:                
                    Domain_file_both.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])) +
                                              '\t' + pair_ortho.split('\t')[1] + '\tboth\n')
                    
                    #check if the considered pair have unmodified domain architecture or one domain lost and store the pair in appropriate file
                    if spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]] == spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]]:
                        Domain_file_nomodif.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])) +
                                                  '\t' + pair_ortho.split('\t')[1] + '\t' + str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])) + '\tnomodif\n')
                    elif len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]]) == len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]]):
                        Domain_file_modif_complex.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])) +
                                                            '\t' + pair_ortho.split('\t')[1] + '\tcomplex_modif\n')
                    else:
                        if len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])+1 == len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]]) or len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])-1 == len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]]):
                            Domain_file_modif_1.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])) +
                                                          '\t' + pair_ortho.split('\t')[1] + '\tmodif\n')
        
                            #check if the loss domain is repeated into the protein or not                    
                            domain_sp1 = spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]]
                            domain_sp2 = spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]]                    
                            
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
                                    
                                    #Inference of the position of the loss
                                    if domain_sp1[1:] == domain_sp2 or domain_sp2[1:] == domain_sp1:
                                        Domain_file_final.write('%s\t%s\t%s\t%s\tf-1\t%s\n' 
                                                                % (pair_ortho.split('\t')[0], str(len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])),
                                                                   pair_ortho.split('\t')[1], str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])),
                                                                   pair_ortho.split('\t')[2]))
                                        
                                    elif domain_sp1[:-1] == domain_sp2 or domain_sp2[:-1] == domain_sp1:
                                        Domain_file_final.write('%s\t%s\t%s\t%s\tb-1\t%s\n' 
                                                                % (pair_ortho.split('\t')[0], str(len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])),
                                                                   pair_ortho.split('\t')[1], str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])),
                                                                   pair_ortho.split('\t')[2]))
                                    
                                    else:
                                        Domain_file_final.write('%s\t%s\t%s\t%s\tint-1\t%s\n' 
                                                                % (pair_ortho.split('\t')[0], str(len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])),
                                                                   pair_ortho.split('\t')[1], str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])),
                                                                   pair_ortho.split('\t')[2]))
                                                                
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
                                    
                                    #Inference of the position of the loss
                                    if domain_sp1[1:] == domain_sp2 or domain_sp2[1:] == domain_sp1:
                                        Domain_file_final.write('%s\t%s\t%s\t%s\tf-1\t%s\n' 
                                                                % (pair_ortho.split('\t')[0], str(len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])),
                                                                   pair_ortho.split('\t')[1], str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])),
                                                                   pair_ortho.split('\t')[2]))
                                        
                                    elif domain_sp1[:-1] == domain_sp2 or domain_sp2[:-1] == domain_sp1:
                                        Domain_file_final.write('%s\t%s\t%s\t%s\tb-1\t%s\n' 
                                                                % (pair_ortho.split('\t')[0], str(len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])),
                                                                   pair_ortho.split('\t')[1], str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])),
                                                                   pair_ortho.split('\t')[2]))
                                    
                                    else:
                                        Domain_file_final.write('%s\t%s\t%s\t%s\tint-1\t%s\n' 
                                                                % (pair_ortho.split('\t')[0], str(len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])),
                                                                   pair_ortho.split('\t')[1], str(len(spec_domain[pair_ortho.split('\t')[1].rstrip('0123456789')][pair_ortho.split('\t')[1]])),
                                                                   pair_ortho.split('\t')[2]))
                                        
                            else:
                                Domain_file_modif_complex.write(pair_ortho)
                                                
                            
                        else:
                            Domain_file_modif_complex.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec_domain[pair_ortho.split('\t')[0].rstrip('0123456789')][pair_ortho.split('\t')[0]])) +
                                                            '\t' + pair_ortho.split('\t')[1] + '\tcomplex_modif\n')
                    
                    
                    
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
    Domain_file_nomodif.close()
    Domain_file_modif_complex.close()
    Domain_file_modif_1.close()
    Domain_file_notrepeated.close()
    Domain_file_repeated.close()
    Domain_file_final.close()
    #
    
#Run the function / with_testis tissue used to infer reference gene
mammals_done = []
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN']):
    for considered_mammals in list_mammals:
        mammals_done.append(considered_mammals)
        for second_mammals in list_mammals:
            if second_mammals not in mammals_done:
                extract_paraortho_modification_group(PairOrtho_in_path = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_paraortho',
                                                     Domain_file_in_path = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/',
                                                     Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paraortho_%s_%s' % (second_mammals, considered_mammals),
                                                     with_testis = True,
                                                     second_species = second_mammals,
                                                     considered_species = considered_mammals)

run_function_in_mammals_dataset()

    
#Run the function / without_testis tissue used to infer reference gene
mammals_done = []
def run_function_in_mammals_dataset(list_mammals = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN']):
    for considered_mammals in list_mammals:
        mammals_done.append(considered_mammals)
        for second_mammals in list_mammals:
            if second_mammals not in mammals_done:
                extract_paraortho_modification_group(PairOrtho_in_path = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_paraortho',
                                                     Domain_file_in_path = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/',
                                                     Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paraortho_%s_%s' % (second_mammals, considered_mammals),
                                                     with_testis = False,
                                                     second_species = second_mammals,
                                                     considered_species = considered_mammals)

run_function_in_mammals_dataset()
