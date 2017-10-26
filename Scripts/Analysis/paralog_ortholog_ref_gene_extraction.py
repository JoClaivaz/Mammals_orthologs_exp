# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:19:02 2017

@author: Claivaz

extract paralog interspecies according to parameter != '1:1' and presence in reference gene.
Reference gene are inferred with Tspec_inference.R for each paralog dataset
only paralog reference gene are assessed in this analysis
"""
def extract_paraortho_refgene_pair_OMA(list_considered_species = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN'],
                                       oma_pair_file = 'D:/UNIL/Master/Master_Project/Data/OMA/oma-pairs.txt',
                                       output_path = 'D:/UNIL/Master/Master_Project/Data/OMA/',
                                       path_reference = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ReferenceGene/',
                                       with_testis_sufix = '_reference_gene_para_notfemale_dataset',
                                       without_testis_sufix = '_reference_gene_para_notfemale_nottestis_dataset',
                                       with_testis = True):
    
    species_done = []
    
    for sp1 in list_considered_species:
        species_done.append(sp1)
        
        for sp2 in list_considered_species:
            if sp2 not in species_done:
                
                if with_testis:
                    ref_gene_sp1 = open('%s%s%s' % (path_reference, sp1, with_testis_sufix), 'r')
                    ref_gene_sp2 = open('%s%s%s' % (path_reference, sp2, with_testis_sufix), 'r')
                
                if not with_testis:
                    ref_gene_sp1 = open('%s%s%s' % (path_reference, sp1, without_testis_sufix), 'r')
                    ref_gene_sp2 = open('%s%s%s' % (path_reference, sp2, without_testis_sufix), 'r')
                    
                ref_gene_list = []
                
                #skip header
                ref_gene_sp1.readline()
                ref_gene_sp2.readline()
                
                for line_tmp in ref_gene_sp1:
                    ref_gene_list.append(line_tmp.split(',')[1].replace('\n', '').replace('"', ''))
                
                for line_tmp in ref_gene_sp2:
                    ref_gene_list.append(line_tmp.split(',')[1].replace('\n', '').replace('"', ''))
                
                ref_gene_sp1.close()
                ref_gene_sp2.close()
                
                if with_testis:
                    paraortho_out = open('%spairwise_paraortho_%s_%s_withtestis' % (output_path, sp1, sp2 ), 'w')
                
                else:
                    paraortho_out = open('%spairwise_paraortho_%s_%s_withouttestis' % (output_path, sp1, sp2 ), 'w')
                
                ortho_file = open(oma_pair_file, 'r')
                
                for ortho_line in ortho_file:
                    
                    if ortho_line.split('\t')[0] in ref_gene_list and ortho_line.split('\t')[1] in ref_gene_list:
                        paraortho_out.write(ortho_line)
                
                ortho_file.close()
                paraortho_out.close()
                
#run function with testis
extract_paraortho_refgene_pair_OMA(with_testis = True)         
#run function without testis
extract_paraortho_refgene_pair_OMA(with_testis = False)