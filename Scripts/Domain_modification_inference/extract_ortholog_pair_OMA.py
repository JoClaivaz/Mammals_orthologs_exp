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
        
        for key_tmp in dict_pair['!1:1'][cons_specie].keys():
            paralog_number += 1
            
            if len(dict_pair['!1:1'][cons_specie][key_tmp]) > 1:
                done_key = []
                list_paralog_tmp = []
                list_paralog_tmp = dict_pair['!1:1'][cons_specie][key_tmp]
                
                for paralog_tmp in list_paralog_tmp:
                    
                    if paralog_tmp not in done_key:
                        done_key.append(paralog_tmp)
                        
                        for paralog_pair in list_paralog_tmp:
                            
                            if paralog_pair not in done_key:
                                
                                try: 
                                    dict_paralog[paralog_tmp]
                                
                                except KeyError: 
                                    dict_paralog[paralog_tmp] = {}
                                
                                try: 
                                    dict_paralog[paralog_tmp][paralog_pair]
                                
                                except KeyError: 
                                    dict_paralog[paralog_tmp][paralog_pair] = [paralog_number, []]
                                    
                                dict_paralog[paralog_tmp][paralog_pair][1].append(key_tmp.rstrip('0123456789'))
           
        for key_tmp in dict_paralog.keys():
            
            for index_tmp in dict_paralog[key_tmp]:
                output_file.write('%s\t%s\tparalog\t%s\t%s\n' % (key_tmp, index_tmp, str(dict_paralog[key_tmp][index_tmp][0]), ' '.join(set(dict_paralog[key_tmp][index_tmp][1]))))
        
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