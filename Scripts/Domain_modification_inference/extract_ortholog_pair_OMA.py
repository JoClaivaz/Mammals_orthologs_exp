# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
171008
This script extract the orthology from OMA file (containing all the pair present in the DB) of the considered species
Create one file ortholog (unique ortholog 1:1) and one file for the paralog of a given specie.

Command line
python extract_ortholog_pair_OMA.py list_considered_species oma_pair_file output_path
"""

def extract_ortholog_pair_OMA(list_considered_species = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN'],
                              oma_pair_file = 'D:/UNIL/Master/Master_Project/Data/OMA/oma-pairs.txt',
                              output_path = 'D:/UNIL/Master/Master_Project/Data/OMA/'):
    
    #store all the pair of interest into dictionnary
    oma_pair = open(oma_pair_file, 'r')
    
    dict_pair = {}
    dict_pair['1:1'] = {}
    dict_pair['!1:1'] = {}
    
    for oma_line in oma_pair:
        
        if oma_line.split('\t')[0].rstrip('0123456789') in list_considered_species and oma_line.split('\t')[1].rstrip('0123456789') in list_considered_species:
            
            if '1:1' in oma_line: 
                
                try:
                    dict_pair['1:1'][oma_line.split('\t')[0].rstrip('0123456789')]
                
                except KeyError:
                    dict_pair['1:1'][oma_line.split('\t')[0].rstrip('0123456789')] = {}
                
                try:
                    dict_pair['1:1'][oma_line.split('\t')[0].rstrip('0123456789')][oma_line.split('\t')[1].rstrip('0123456789')]
                
                except KeyError:
                    dict_pair['1:1'][oma_line.split('\t')[0].rstrip('0123456789')][oma_line.split('\t')[1].rstrip('0123456789')] = []
                    
                    
                dict_pair['1:1'][oma_line.split('\t')[0].rstrip('0123456789')][oma_line.split('\t')[1].rstrip('0123456789')].append(oma_line)
            
            else:
                
                try:
                    dict_pair['!1:1'][oma_line.split('\t')[0].rstrip('0123456789')]
                
                except KeyError:
                    dict_pair['!1:1'][oma_line.split('\t')[0].rstrip('0123456789')] = {}    
                
                try:
                    dict_pair['!1:1'][oma_line.split('\t')[1].rstrip('0123456789')]
                
                except KeyError:
                    dict_pair['!1:1'][oma_line.split('\t')[1].rstrip('0123456789')] = {}
                
                try:
                    dict_pair['!1:1'][oma_line.split('\t')[0].rstrip('0123456789')][oma_line.split('\t')[1]]
                
                except KeyError:
                    dict_pair['!1:1'][oma_line.split('\t')[0].rstrip('0123456789')][oma_line.split('\t')[1]] = []
                
                try:
                    dict_pair['!1:1'][oma_line.split('\t')[1].rstrip('0123456789')][oma_line.split('\t')[0]]
                
                except KeyError:
                    dict_pair['!1:1'][oma_line.split('\t')[1].rstrip('0123456789')][oma_line.split('\t')[0]] = []
                    
                dict_pair['!1:1'][oma_line.split('\t')[1].rstrip('0123456789')][oma_line.split('\t')[0]].append(oma_line.split('\t')[1])
                dict_pair['!1:1'][oma_line.split('\t')[0].rstrip('0123456789')][oma_line.split('\t')[1]].append(oma_line.split('\t')[0])    
                    
    oma_pair.close() 
    #
    
    #generation of paralog pair into one specie file
    for cons_specie in list_considered_species:
        output_file = open('%spairwise_paralog_%s' % (output_path, cons_specie), 'w')
        dict_paralog = {}
        paralog_number = 0
        fused_family = {}
        
        for key_tmp in dict_pair['!1:1'][cons_specie].keys():
            paralog_number += 1
            
            if len(dict_pair['!1:1'][cons_specie][key_tmp]) > 1:
                list_paralog_tmp = []
                list_paralog_tmp = dict_pair['!1:1'][cons_specie][key_tmp]
                any_in_dict = False
                considered_family = 0
                
                for paralog_tmp in list_paralog_tmp:
                    try:
                        dict_paralog[paralog_tmp]
                        any_in_dict = True
                        
                        if considered_family == 0:
                            considered_family = dict_paralog[paralog_tmp]
                        elif considered_family != dict_paralog[paralog_tmp]:
                            try:
                                fused_family[considered_family]
                                fused_family[considered_family].append(dict_paralog[paralog_tmp])
                            except KeyError:
                                fused_family[considered_family] = [dict_paralog[paralog_tmp]]
                            try:
                                fused_family[dict_paralog[paralog_tmp]]
                                fused_family[dict_paralog[paralog_tmp]].append(considered_family)
                            except KeyError:
                                fused_family[dict_paralog[paralog_tmp]] = [considered_family]
                    except KeyError:
                        pass
                
                for paralog_tmp in list_paralog_tmp:
                    if any_in_dict == False:
                        dict_paralog[paralog_tmp] = paralog_number
                        
                    else:
                        try:
                            dict_paralog[paralog_tmp]
                        except KeyError:
                            dict_paralog[paralog_tmp] = considered_family
                        
        
        dict_family = {}            
        for key_tmp in dict_paralog.keys():
            try:
                dict_family[dict_paralog[key_tmp]].append(key_tmp)
            except KeyError:
                dict_family[dict_paralog[key_tmp]] = [key_tmp]
                
        done_family = []
        for key_tmp in fused_family.keys():
            fused_numbers = fused_family[key_tmp]
            
            for cons_number in fused_numbers:
                if cons_number not in done_family and len(dict_family[cons_number]) > 0:
                    
                    for new_paralog in dict_family[cons_number]:
                        if new_paralog not in dict_family[key_tmp]:
                            dict_family[key_tmp].append(new_paralog)
                        
                    done_family.append(cons_number)
                    dict_family[cons_number] = []
            
                    
           
        for key_tmp in dict_family.keys():
            gene_done = []
            
            for gene_tmp in dict_family[key_tmp]:
                gene_done.append(gene_tmp)
                for gene_2 in dict_family[key_tmp]:
                    if gene_2 not in gene_done:
                        output_file.write('%s\t%s\tparalog\t%s\t\n' % (gene_tmp, gene_2, str(key_tmp)))
        
        output_file.close()
    #
    
    #generation of unique ortholog pair
    specie_done = []
    for cons_specie in list_considered_species:
        specie_done.append(cons_specie)
        for specie_tmp in list_considered_species:
            if specie_tmp not in specie_done:
                output_file = open('%spairwise_ortholog_%s_%s' % (output_path, cons_specie, specie_tmp), 'w')
                try:
                    dict_pair['1:1'][cons_specie][specie_tmp]
                    output_file.write(''.join(dict_pair['1:1'][cons_specie][specie_tmp]))
                except KeyError:
                    pass
                try:
                    dict_pair['1:1'][specie_tmp][cons_specie]
                    output_file.write(''.join(dict_pair['1:1'][specie_tmp][cons_specie]))
                except KeyError:
                    pass
                
                output_file.close()
    #
                     
#Run FUN in command line
#python extract_ortholog_pair_OMA.py list_considered_species oma_pair_file output_path
import sys
extract_ortholog_pair_OMA(list_considered_species = sys.argv[1],
                          oma_pair_file = sys.argv[2],
                          output_path = sys.argv[3])                    

#In python
#extract_ortholog_pair_OMA()