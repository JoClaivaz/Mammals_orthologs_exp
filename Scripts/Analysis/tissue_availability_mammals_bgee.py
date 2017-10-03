# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 17:12:07 2017

@author: Claivaz
Extract the different considered tissue available in Bgee RNAseq files, allowing the selection of the different considered species
"""
import os

def extract_unique_tissue(species_dir, output_file):
    different_species = os.listdir(species_dir)
    dict_species = {}
    
    for considered_species in different_species:
        species_tissue = []
        different_dataset = os.listdir(species_dir + considered_species)
        
        for considered_file in different_dataset:
            if '.zip' not in considered_file:
                open_file = open('%s%s/%s' % (species_dir, considered_species, considered_file), 'r')
                open_file.readline()
                
                for line in open_file:
                    if (line.split('\t')[5].replace('"', '') + '_' + line.split('\t')[7].replace('"', '')) not in species_tissue:
                        species_tissue.append('%s_%s' % (line.split('\t')[5].replace('"', ''), line.split('\t')[7].replace('"', '')))
                
                open_file.close()
        dict_species[considered_species] = species_tissue            
    
    output_open = open(output_file, 'w')
    
    for considered_key in dict_species.keys():
        output_open.write('%s\n%s\n\n' % (considered_key, '\n'.join(dict_species[considered_key])))
    
    output_open.close()
        
extract_unique_tissue(species_dir = 'E:/Bgee/', 
                      output_file = 'D:/UNIL/Master/Master_Project/Data/Bgee/tissue_available_species_bgee')

