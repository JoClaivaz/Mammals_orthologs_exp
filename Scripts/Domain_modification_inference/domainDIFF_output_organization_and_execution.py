# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 20:28:30 2017

@author: jclaivaz

this script create the pair domain file mendatory for domainDIFF, and execute the domainDIFF program with the output
"""

def domainDIFF_output_organization_and_execution(list_species = ['BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX'], central_species = 'HUMAN', path_input= '/media/jclaivaz/Data/UNIL/Master/Master_Project/Data/domain_architecture_inference/', regexp_output = 'ortholog', path_domainDIFF = '/home/jclaivaz/Bureau/PfamScan/domainDiff'):
	
	#organization
	central_species_file = open('%s%s_domain' % (path_input, central_species), 'r')
	central_species_text = ''

	for central_species_line in central_species_file:
		if '#' not in central_species_line:
			central_species_text += central_species_line

	central_species_file.close()

	for considered_species in list_species:
		considered_file = open('%s%s_domain' % (path_input, considered_species), 'r')
		output_file_considered = '%s%s_%s_%s_domain' % (path_input, regexp_output, central_species, considered_species)
		output_file = open(output_file_considered, 'w')

		for considered_line in considered_file:
			output_file.write(considered_line)

		output_file.write(central_species_text.replace('\nHUMAN00008', 'HUMAN00008'))

		output_file.close()
		considered_file.close()
		#

		#exececute domainDIFF
		import subprocess
		
		bash_command = '%s -a %s > %sortholog_%s_%s_domain_modifications' % (path_domainDIFF, output_file_considered, path_input, central_species, considered_species)
		print(bash_command)
	
		#DOESN'T WORK REWORK ON IT
		#subprocess.run(bash_command.split(' '), shell = True)
		#

###RUN function
#on linux
domainDIFF_output_organization_and_execution()