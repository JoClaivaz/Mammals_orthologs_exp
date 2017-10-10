# -*- coding: utf-8 -*-
"""
171002, last modification 171010

@author: Claivaz
parsed the expression files: 
    -keeping only comparable state
    -Gene present in 1 domain not repeated or no modification
        -store all the pairwise expression data 
    -use ENSEMBL identifier and replace it with the OMA one
    -create one expression file for each mammals species

N.B.: expression directories is organized as one folder per species containing the specific expression files 
and domain inference files are organized as it was expected through the domain inference pipeline

Command line
python extraction_state_expression_file.py considered_tissues considered_species control_modification_pairs_folder_path  path_expression_dir oma_ensembl_converter output_file_path
"""

def extraction_state_expression_file(considered_tissues = ['adult mammalian kidney', 'testis', 'heart', 'brain', 'skeletal muscle tissue',
                                                           'colon', 'lung', 'spleen', 'cerebellum', 'liver', 'frontal cortex',
                                                           'prefrontal cortex', 'female gonad', 'placenta', 'kidney'],
                                     considered_species = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN'],
                                     control_modification_pairs_folder_path = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/',
                                     path_expression_dir = 'E:/Bgee/',
                                     oma_ensembl_converter = 'D:/UNIL/Master/Master_Project/Data/OMA/fasta_seq/oma-ensembl.txt',
                                     output_file_path = 'D:/UNIL/Master/Master_Project/Data/Bgee/'):
                                     #prefix_ensembl_ID = ['ENSBTG', 'ENSGG', 'ENSG', 'ENSMMUG', 'ENSMODG', 'ENSMUSG', 'ENSPTRG', 'ENSSSCG', 'ENSRNOG'] not necessary

    #store considered gene (1 domain modif / control) into list for each species
    dict_gene = {}
    for considered_species_tmp in considered_species:
        dict_gene[considered_species_tmp] = {}
    
    for central_species in considered_species: 
        for considered_species_tmp in considered_species:
            
            control_modif_pair_name = ['ortholog_%s_%s_domain_nomodif' % (central_species, considered_species_tmp),
                                       'putative_ortholog_%s_%s_domain_loss' % (central_species, considered_species_tmp)]
            
            for control_or_modif in control_modif_pair_name:
                
                try:
                    open_file_conmodif = open(control_modification_pairs_folder_path + control_or_modif, 'r')
                    
                    for line_conmodif in open_file_conmodif:
                        
                        dict_gene[line_conmodif.split('\t')[0].rstrip('0123456789')][line_conmodif.split('\t')[0]] = ''
                        dict_gene[line_conmodif.split('\t')[2].rstrip('0123456789')][line_conmodif.split('\t')[2]] = ''
                        
                    open_file_conmodif.close()
                    
                except:
                    pass
        
        control_modif_pair_name = ['%sputative_paralog_%s_domain_loss' % (control_modification_pairs_folder_path, central_species),
                                   '%sparalog_%s_domain_nomodif' % (control_modification_pairs_folder_path, central_species)] 
        
        for conmod_file in control_modif_pair_name:
            open_file_conmodif = open(conmod_file, 'r')
            
            for line_conmodif in open_file_conmodif:
                        
                dict_gene[line_conmodif.split('\t')[0].rstrip('0123456789')][line_conmodif.split('\t')[0]] = ''
                dict_gene[line_conmodif.split('\t')[2].rstrip('0123456789')][line_conmodif.split('\t')[2]] = ''
            
            open_file_conmodif.close()
        
    #store ensembl and oma identifier considered into dictionary
    omaensembl_file = open(oma_ensembl_converter, 'r')
    dict_converter = {}
    
    for omaensembl_line in omaensembl_file:
        oma_identifier = omaensembl_line.split('\t')[0] 
        
        try:
            dict_gene[oma_identifier.rstrip('0123456789')][oma_identifier]
            dict_converter[omaensembl_line.split('\t')[1].split('.')[0].split('\n')[0]] = oma_identifier
            
        except KeyError:
            pass
        
    omaensembl_file.close()
    
    #parsing expression file in function of the considered pair gene and tissues
    import os
    for all_species in considered_species:
        
        different_dataset = os.listdir(path_expression_dir + all_species)
        expression_output = open('%s%s_expression_parsed' % (output_file_path, all_species), 'w')
        expression_output.write('ExperimentID\tLibraryID\tGeneID\tAnatomicalEntityName\tStageName\tSex\tReadCount\tTPM\tFPKM\n')
        
        for considered_file in different_dataset:
            if '.zip' not in considered_file:
                open_file = open('%s%s/%s' % (path_expression_dir, all_species, considered_file), 'r')
     
                for line_exp in open_file:
                    try:
                        dict_converter[line_exp.split('\t')[3].split('.')[0]]
                        
                        if line_exp.split('\t')[5].replace('"', '') in considered_tissues:
                            expression_output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' 
                                                    % (line_exp.split('\t')[0], line_exp.split('\t')[1], 
                                                       dict_converter[line_exp.split('\t')[3].split('.')[0]],
                                                       line_exp.split('\t')[5].replace('"', ''), line_exp.split('\t')[7], 
                                                       line_exp.split('\t')[8], line_exp.split('\t')[10],
                                                       line_exp.split('\t')[11], line_exp.split('\t')[12]))
                        
                    except KeyError:
                        pass
        
                open_file.close()
        expression_output.close()
            
#Run FUN in command line
#python extraction_state_expression_file.py considered_tissues considered_species control_modification_pairs_folder_path  path_expression_dir oma_ensembl_converter output_file_path
import sys
extraction_state_expression_file(considered_tissues = sys.argv[1],
                                 considered_species = sys.argv[2],
                                 control_modification_pairs_folder_path = sys.argv[3],
                                 path_expression_dir = sys.argv[4],
                                 oma_ensembl_converter = sys.argv[5],
                                 output_file_path = sys.argv[6])
#Run in python
#extraction_state_expression_file()