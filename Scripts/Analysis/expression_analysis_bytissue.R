'''
Joaquim Claivaz
171023

Expression analysis plot results
see difference in evolution of the conservation expression in function of domain modification
'''
####FUN####
data_organization_ExpAnal_ortholog = function(considered_species_name_1, 
                                            considered_species_name_2,
                                            considered_sex_vector = F,
                                            notconsidered_sex_vector = c('female'),
                                            notconsidered_anat_vector = c('kidney', 'frontal cortex'),
                                            considered_anat_vector = F,
                                            considered_devtime_vector = F,
                                            notconsidered_devtime_vector = c("9th week post-fertilization human stage (human)",
                                                                            "10th week post-fertilization human stage (human)",
                                                                            "16th week post-fertilization human stage (human)",
                                                                            "17th week post-fertilization human stage (human)",
                                                                            "19th week post-fertilization human stage (human)",
                                                                            "9th week post-fertilization human stage (human)",
                                                                            "16th week post-fertilization human stage (human)",
                                                                            "19th week post-fertilization human stage (human)",
                                                                            "9th week post-fertilization human stage (human)",
                                                                            "15th week post-fertilization human stage (human)",
                                                                            "19th week post-fertilization human stage (human)"),
                                            expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                            domain_modif_sufix = '_domain_loss'){
  
  #open_files
  specie1_data = read.table(paste0(expression_data_path_prefix,
                                   considered_species_name_1,
                                   expression_data_sufix), header = T, sep = '\t')
  specie2_data = read.table(paste0(expression_data_path_prefix,
                                   considered_species_name_2,
                                   expression_data_sufix), header = T, sep = '\t')
  
  options(show.error.messages = FALSE)
  if (class(try(read.table(paste0(domain_control_path_prefix,
                                  considered_species_name_2, '_', considered_species_name_1,
                                  domain_control_sufix), header = F, sep = '\t'))) != 'try-error'){
    domain_control = read.table(paste0(domain_control_path_prefix,
                                       considered_species_name_2, '_', considered_species_name_1,
                                       domain_control_sufix), header = F, sep = '\t')
    domain_modif = read.table(paste0(domain_modif_path_prefix,
                                     considered_species_name_2, '_', considered_species_name_1,
                                     domain_modif_sufix), header = F, sep = '\t')
  }else{
    domain_control = read.table(paste0(domain_control_path_prefix,
                                       considered_species_name_1, '_', considered_species_name_2,
                                       domain_control_sufix), header = F, sep = '\t')
    domain_modif = read.table(paste0(domain_modif_path_prefix,
                                     considered_species_name_1, '_', considered_species_name_2,
                                     domain_modif_sufix), header = F, sep = '\t')
  }
  options(show.error.messages = TRUE)
  #
  
  #Sex
  if (considered_sex_vector != F){
    specie1_data = specie1_data[specie1_data$Sex %in% considered_sex_vector,]
    specie2_data = specie2_data[specie2_data$Sex %in% considered_sex_vector,]
  } 
  if (notconsidered_sex_vector != F){
    specie1_data = specie1_data[!(specie1_data$Sex %in% notconsidered_sex_vector),]
    specie2_data = specie2_data[!(specie2_data$Sex %in% notconsidered_sex_vector),]
    
  }
  if (notconsidered_sex_vector == F & considered_sex_vector == F){
    notkeep_Sex = unique(specie1_data$Sex)[!(unique(specie1_data$Sex) %in% unique(specie2_data$Sex))]
    specie1_data = specie1_data[!(specie1_data$Sex %in% notkeep_Sex),]
    specie2_data = specie2_data[!(specie2_data$Sex %in% notkeep_Sex),]
  }
  #
  
  #AnatomicalEntityName
  keep_AnatomicalEntityName = unique(specie2_data$AnatomicalEntityName)[unique(specie2_data$AnatomicalEntityName) %in% unique(specie1_data$AnatomicalEntityName)]
  specie1_data = specie1_data[specie1_data$AnatomicalEntityName %in% keep_AnatomicalEntityName,]
  specie2_data = specie2_data[specie2_data$AnatomicalEntityName %in% keep_AnatomicalEntityName,]
  # keep_AnatomicalEntityName = unique(specie1_data$AnatomicalEntityName)[unique(specie1_data$AnatomicalEntityName) %in% unique(specie2_data$AnatomicalEntityName)]
  # specie1_data = specie1_data[specie1_data$AnatomicalEntityName %in% keep_AnatomicalEntityName,]
  # specie2_data = specie2_data[specie2_data$AnatomicalEntityName %in% keep_AnatomicalEntityName,]
  #
  if (considered_anat_vector != F){
    specie1_data = specie1_data[specie1_data$AnatomicalEntityName %in% considered_anat_vector,]
    specie2_data = specie2_data[specie2_data$AnatomicalEntityName %in% considered_anat_vector,]
    
  }
  if (notconsidered_anat_vector != F){
    specie1_data = specie1_data[!(specie1_data$AnatomicalEntityName %in% notconsidered_anat_vector),]
    specie2_data = specie2_data[!(specie2_data$AnatomicalEntityName %in% notconsidered_anat_vector),]
  }
  
  #developmental time filtering
  if (considered_devtime_vector != F){
    specie1_data = specie1_data[specie1_data$StageName %in% considered_devtime_vector,]
    specie2_data = specie2_data[specie2_data$StageName %in% considered_devtime_vector,]
    
  }
  if (notconsidered_devtime_vector != F){
    specie1_data = specie1_data[!(specie1_data$StageName %in% notconsidered_devtime_vector),]
    specie2_data = specie2_data[!(specie2_data$StageName %in% notconsidered_devtime_vector),]
  }
  #
  
  #MEAN tpm in function of geneID and AnatomicalEntityName / Sex not considered (for each sex the function should be run separatedly)
  specie2_data = aggregate(specie2_data$TPM ~ specie2_data$GeneID + as.character(specie2_data$AnatomicalEntityName), FUN = mean)
  names(specie2_data) = c('GeneID', 'Tissue', 'TPM')
  specie1_data = aggregate(specie1_data$TPM ~ specie1_data$GeneID + as.character(specie1_data$AnatomicalEntityName), FUN = mean)
  names(specie1_data) = c('GeneID', 'Tissue', 'TPM')
  #
  
  #group the domain inference data
  domain_modif$V6 = NULL
  domain_status = rbind(domain_control, domain_modif)
  names(domain_status) = c('GeneID_1', 'len_1', 'GeneID_2', 'len_2', 'status')
  domain_status$status = as.character(domain_status$status)
  domain_status$status[domain_status$status == 'nomodif'] = 'control'
  #
  
  #convert tpm into numerics and log transformation
  specie1_data$TPM = log2(as.numeric(specie1_data$TPM) + 0.000001)
  specie2_data$TPM = log2(as.numeric(specie2_data$TPM) + 0.000001)
  specie1_data$TPM[specie1_data$TPM < 1] = 0  
  specie2_data$TPM[specie2_data$TPM < 1] = 0
  #
  
  #domain_status file numerotation of pair
  domain_status$pair_nbr = 1 : dim(domain_status)[1]
  #
  
  #add pair number and status in species files
  specie1_data$pair_nbr = NA
  specie1_data$status = NA
  for (row_tmp in 1:dim(specie1_data)[1]){
    if(specie1_data$GeneID[row_tmp] %in% domain_status$GeneID_1 | specie1_data$GeneID[row_tmp] %in% domain_status$GeneID_2){
      specie1_data$pair_nbr[row_tmp] = domain_status$pair_nbr[domain_status$GeneID_1 == as.character(specie1_data$GeneID[row_tmp]) | domain_status$GeneID_2 == as.character(specie1_data$GeneID[row_tmp])]
      specie1_data$status[row_tmp] = domain_status$status[domain_status$GeneID_1 == as.character(specie1_data$GeneID[row_tmp]) | domain_status$GeneID_2 == as.character(specie1_data$GeneID[row_tmp])]
    }
  }
  specie2_data$pair_nbr = NA
  specie2_data$status = NA
  for (row_tmp in 1:dim(specie2_data)[1]){
    if(specie2_data$GeneID[row_tmp] %in% domain_status$GeneID_1 | specie2_data$GeneID[row_tmp] %in% domain_status$GeneID_2){
      specie2_data$pair_nbr[row_tmp] = domain_status$pair_nbr[domain_status$GeneID_1 == as.character(specie2_data$GeneID[row_tmp]) | domain_status$GeneID_2 == as.character(specie2_data$GeneID[row_tmp])]
      specie2_data$status[row_tmp] = domain_status$status[domain_status$GeneID_1 == as.character(specie2_data$GeneID[row_tmp]) | domain_status$GeneID_2 == as.character(specie2_data$GeneID[row_tmp])]
    }
  }
  #
  
  #keep complete case
  specie1_data = specie1_data[complete.cases(specie1_data),]
  specie2_data = specie2_data[complete.cases(specie2_data),]
  #
  
  #keep_name of each dataset
  specie1_name = gsub('[[:digit:]]', '', specie1_data$GeneID[1])
  specie2_name = gsub('[[:digit:]]', '', specie2_data$GeneID[1])
  #
  
  #check if pair are presented in both dataset
  keep_pair = unique(specie2_data$pair_nbr)[unique(specie2_data$pair_nbr) %in% unique(specie1_data$pair_nbr)]
  specie1_data = specie1_data[specie1_data$pair_nbr %in% keep_pair,]
  specie2_data = specie2_data[specie2_data$pair_nbr %in% keep_pair,]
  #
  
  #group data
  names_tissue = unique(as.character(specie2_data$Tissue))
  specie1_data = aggregate(TPM ~ pair_nbr + status, data = specie1_data, paste, collapse = '_')
  specie1_data = separate(specie1_data, col = 'TPM', sep= '_', into = names_tissue)
  specie2_data = aggregate(TPM ~ pair_nbr + status, data = specie2_data, paste, collapse = '_')
  specie2_data = separate(specie2_data, col = 'TPM', sep= '_', into = names_tissue)
  #
  
  return(list(specie1_data, specie1_name, specie2_data, specie2_name))
}

data_organization_ExpAnal_paralog = function(considered_species_name,
                                             considered_sex_vector = F,
                                             notconsidered_sex_vector = c('female'),
                                             notconsidered_anat_vector = c('kidney', 'frontal cortex'),
                                             considered_anat_vector = F,
                                             considered_devtime_vector = F,
                                             notconsidered_devtime_vector = c("9th week post-fertilization human stage (human)",
                                                                              "10th week post-fertilization human stage (human)",
                                                                              "16th week post-fertilization human stage (human)",
                                                                              "17th week post-fertilization human stage (human)",
                                                                              "19th week post-fertilization human stage (human)",
                                                                              "9th week post-fertilization human stage (human)",
                                                                              "16th week post-fertilization human stage (human)",
                                                                              "19th week post-fertilization human stage (human)",
                                                                              "9th week post-fertilization human stage (human)",
                                                                              "15th week post-fertilization human stage (human)",
                                                                              "19th week post-fertilization human stage (human)"),
                                             expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                             expression_data_sufix = '_expression_parsed',
                                             domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                             domain_control_sufix = '_domain_nomodif',
                                             domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                             domain_modif_sufix = '_domain_loss'){
  
  #open_files
  specie_data = read.table(paste0(expression_data_path_prefix,
                                  considered_species_name,
                                  expression_data_sufix), header = T, sep = '\t')
  
  domain_control = read.table(paste0(domain_control_path_prefix,
                                     considered_species_name,
                                     domain_control_sufix), header = F, sep = '\t')
  domain_modif = read.table(paste0(domain_modif_path_prefix,
                                   considered_species_name,
                                   domain_modif_sufix), header = F, sep = '\t')
  #
  
  #Sex
  if (considered_sex_vector != F){
    specie_data = specie_data[specie_data$Sex %in% considered_sex_vector,]
  } 
  if (notconsidered_sex_vector != F){
    specie_data = specie_data[!(specie_data$Sex %in% notconsidered_sex_vector),]
  }
  #
  
  #AnatomicalEntityName
  if (considered_anat_vector != F){
    specie_data = specie_data[specie_data$AnatomicalEntityName %in% considered_anat_vector,]
  }
  if (notconsidered_anat_vector != F){
    specie_data = specie_data[!(specie_data$AnatomicalEntityName %in% notconsidered_anat_vector),]
  }
  
  #developmental time filtering
  if (considered_devtime_vector != F){
    specie_data = specie_data[specie_data$StageName %in% considered_devtime_vector,]
  }
  if (notconsidered_devtime_vector != F){
    specie_data = specie_data[!(specie_data$StageName %in% notconsidered_devtime_vector),]
  }
  #
  
  #MEAN tpm in function of geneID and AnatomicalEntityName / Sex not considered (for each sex the function should be run separatedly)
  specie_data = aggregate(specie_data$TPM ~ specie_data$GeneID + specie_data$AnatomicalEntityName, FUN = mean)
  names(specie_data) = c('GeneID', 'Tissue', 'TPM')
  #
  
  #convert tpm into numerics and log transformation
  specie_data$TPM = log2(as.numeric(specie_data$TPM) + 0.000001)
  specie_data$TPM[specie_data$TPM < 1] = 0
  #
  
  #group the domain inference data
  domain_control$V8 = NULL
  domain_status = rbind(domain_control, domain_modif)
  domain_status$V5 = as.character(domain_status$V5)
  domain_status$V5[domain_status$V5 == 'paralog'] = 'control'
  names(domain_status) = c('GeneID_1', 'len_1', 'GeneID_2', 'len_2', 'status', 'ParalogGroup', 'SpeciesParalogInference')
  domain_status = aggregate(ParalogGroup ~ GeneID_1 + GeneID_2 + len_1 + len_2, data = domain_status, FUN = paste0)
  #
  
  return(list(specie_data, domain_status))
}
####

require(tidyr)

####HUMAN BOVIN orthologs####
{
list_output = data_organization_ExpAnal_ortholog(considered_species_name_1 = 'BOVIN', 
                                                 considered_species_name_2 = 'HUMAN')
specie1_data = list_output[[1]]
specie1_name = list_output[[2]]
specie2_data = list_output[[3]]
specie2_name = list_output[[4]]

#test if the order of pair number are the same
all(specie1_data$pair_nbr == specie2_data$pair_nbr)
#

#plot
tissue_list = colnames(specie1_data)[3:length(colnames(specie1_data))]

{
pdf(paste0('C:/Users/Claivaz/Desktop/', 'expression_analysis_HUMAN_BOVIN' ,'.pdf'))
par(mfrow = c(5,4), mai=c(0.4,0.35,0.3,0.01))
  
for (tissue_tmp in 1:length(tissue_list)){
  
  smoothScatter(specie1_data[,(tissue_tmp + 2)][specie1_data$status == 'control'], 
                specie2_data[,(tissue_tmp + 2)][specie2_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = paste0(tissue_list[tissue_tmp], ' - control'), cex.main = 0.8, cex.lab = 0.6, xlim = c(0, 20), ylim = c(0, 20))
  title(xlab = paste0(specie1_name, ' expression level (log transformed)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(specie2_name, ' expression level (log transformed)'), cex.lab =0.6, line = 2)
  value_tmp1 = as.numeric(specie1_data[,(tissue_tmp + 2)][specie1_data$status == 'control'])
  value_tmp2 = as.numeric(specie2_data[,(tissue_tmp + 2)][specie2_data$status == 'control'])
  abline(lm(value_tmp2 ~ value_tmp1), col = 'red')
  linear_param = lm(value_tmp2 ~ value_tmp1)$coefficients
  cor_value = cor.test(value_tmp2, value_tmp1)$estimate[[1]]
  text(10, 5, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(specie1_data[,(tissue_tmp + 2)][specie1_data$status != 'control'], 
                specie2_data[,(tissue_tmp + 2)][specie2_data$status != 'control'],
                xlab = '',
                ylab = '',
                main = paste0(tissue_list[tissue_tmp], ' - domain modification'), cex.main = 0.8, cex.lab = 0.6, xlim = c(0, 20), ylim = c(0, 20))
  title(xlab = paste0(specie1_name, ' expression level (log transformed)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(specie2_name, ' expression level (log transformed)'), cex.lab =0.6, line = 2)
  value_tmp1 = as.numeric(specie1_data[,(tissue_tmp + 2)][specie1_data$status != 'control'])
  value_tmp2 = as.numeric(specie2_data[,(tissue_tmp + 2)][specie2_data$status != 'control'])
  abline(lm(value_tmp2 ~ value_tmp1), col = 'red')
  linear_param = lm(value_tmp2 ~ value_tmp1)$coefficients
  cor_value = cor.test(value_tmp2, value_tmp1)$estimate[[1]]
  text(10, 5, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
}

dev.off()
}
}
#

####HUMAN MOUSE orthologs####
{
list_output = data_organization_ExpAnal_ortholog(considered_species_name_1 = 'MOUSE', 
                                                 considered_species_name_2 = 'HUMAN')
specie1_data = list_output[[1]]
specie1_name = list_output[[2]]
specie2_data = list_output[[3]]
specie2_name = list_output[[4]]

#test if the order of pair number are the same
all(specie1_data$pair_nbr == specie2_data$pair_nbr)
#

#plot
tissue_list = colnames(specie1_data)[3:length(colnames(specie1_data))]

{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'expression_analysis_HUMAN_MOUSE' ,'.pdf'))
  par(mfrow = c(5,4), mai=c(0.4,0.35,0.3,0.01))
  
  for (tissue_tmp in 1:length(tissue_list)){
    
    smoothScatter(specie1_data[,(tissue_tmp + 2)][specie1_data$status == 'control'], 
                  specie2_data[,(tissue_tmp + 2)][specie2_data$status == 'control'],
                  xlab = '',
                  ylab = '',
                  main = paste0(tissue_list[tissue_tmp], ' - control'), cex.main = 0.8, cex.lab = 0.6, xlim = c(0, 20), ylim = c(0, 20))
    title(xlab = paste0(specie1_name, ' expression level (log transformed)'), cex.lab =0.6, line = 2)
    title(ylab = paste0(specie2_name, ' expression level (log transformed)'), cex.lab =0.6, line = 2)
    value_tmp1 = as.numeric(specie1_data[,(tissue_tmp + 2)][specie1_data$status == 'control'])
    value_tmp2 = as.numeric(specie2_data[,(tissue_tmp + 2)][specie2_data$status == 'control'])
    abline(lm(value_tmp2 ~ value_tmp1), col = 'red')
    linear_param = lm(value_tmp2 ~ value_tmp1)$coefficients
    cor_value = cor.test(value_tmp2, value_tmp1)$estimate[[1]]
    text(10, 5, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
    
    smoothScatter(specie1_data[,(tissue_tmp + 2)][specie1_data$status != 'control'], 
                  specie2_data[,(tissue_tmp + 2)][specie2_data$status != 'control'],
                  xlab = '',
                  ylab = '',
                  main = paste0(tissue_list[tissue_tmp], ' - domain modification'), cex.main = 0.8, cex.lab = 0.6, xlim = c(0, 20), ylim = c(0, 20))
    title(xlab = paste0(specie1_name, ' expression level (log transformed)'), cex.lab =0.6, line = 2)
    title(ylab = paste0(specie2_name, ' expression level (log transformed)'), cex.lab =0.6, line = 2)
    value_tmp1 = as.numeric(specie1_data[,(tissue_tmp + 2)][specie1_data$status != 'control'])
    value_tmp2 = as.numeric(specie2_data[,(tissue_tmp + 2)][specie2_data$status != 'control'])
    abline(lm(value_tmp2 ~ value_tmp1), col = 'red')
    linear_param = lm(value_tmp2 ~ value_tmp1)$coefficients
    cor_value = cor.test(value_tmp2, value_tmp1)$estimate[[1]]
    text(10, 5, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
    
  }
  
  dev.off()
}
}
#

####HUMAN paralogs####
{
list_output = data_organization_ExpAnal_paralog(considered_species_name = 'HUMAN')
specie_data = list_output[[1]]
domain_status = list_output[[2]]
{
pdf(paste0('C:/Users/Claivaz/Desktop/', 'expression_analysis_HUMAN' ,'.pdf'))
par(mfrow = c(6,3), mai=c(0.4,0.35,0.3,0.01))
  
tissue_list = as.character(unique(specie_data$Tissue))
paralog_family = unique(domain_status$ParalogGroup)

#Infer the group status and defined as gene1 the longer pair in modification group
domain_status$GeneID_1 = as.character(domain_status$GeneID_1)
domain_status$GeneID_2 = as.character(domain_status$GeneID_2)
domain_status$status = 'control'
domain_status$status[domain_status$len_1 != domain_status$len_2] = 'modif'

for(row_tmp in 1:dim(domain_status)[1]){
  if(domain_status$len_1[row_tmp] < domain_status$len_2[row_tmp]){
    tmp_GeneID = domain_status$GeneID_1[row_tmp]
    domain_status$GeneID_1[row_tmp] = domain_status$GeneID_2[row_tmp]
    domain_status$GeneID_2[row_tmp] = tmp_GeneID
  }
}
#

for (tissue_tmp in 1:length(tissue_list)){
  specie_tmp = specie_data[specie_data$Tissue == tissue_list[tissue_tmp],]
  domain_tmp = domain_status
  
  specie_tmp$ParalogGroup = NA
  for (row_tmp in 1:dim(specie_tmp)[1]){
    specie_tmp$ParalogGroup[row_tmp] = as.numeric(domain_status$ParalogGroup[domain_status$GeneID_1 == as.character(specie_tmp$GeneID[row_tmp]) | domain_status$GeneID_2 == as.character(specie_tmp$GeneID[row_tmp])][1])
  }
  
  specie_tmp = specie_tmp[complete.cases(specie_tmp),]
  
  #generate reference gene per paralog family
  reference_gene = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/ReferenceGene/', gsub('[[:digit:]]', '', domain_control[1,1])))
  names(reference_gene)[2] = 'GeneID'
  
  domain_tmp$exp_1 = NA
  domain_tmp$exp_2 = NA
  domain_modif = domain_tmp[domain_tmp$status == 'modif',]
  domain_control = domain_tmp[domain_tmp$status == 'control',]
  
  #only maximal expression gene considered as reference gene and are kept in geneID1
  for(row_tmp in 1:dim(domain_control)[1]){
    if(domain_control$GeneID_1[row_tmp] %in% reference_gene$GeneID & as.character(domain_control$GeneID_2[row_tmp]) %in% specie_tmp$GeneID){
      domain_control$exp_1[row_tmp] = specie_tmp$TPM[specie_tmp$GeneID == as.character(domain_control$GeneID_1[row_tmp])]
      domain_control$exp_2[row_tmp] = specie_tmp$TPM[specie_tmp$GeneID == as.character(domain_control$GeneID_2[row_tmp])]
    }else if(domain_control$GeneID_2[row_tmp] %in% reference_gene$GeneID & as.character(domain_control$GeneID_1[row_tmp]) %in% specie_tmp$GeneID){
      domain_control$exp_2[row_tmp] = specie_tmp$TPM[specie_tmp$GeneID == as.character(domain_control$GeneID_1[row_tmp])]
      domain_control$exp_1[row_tmp] = specie_tmp$TPM[specie_tmp$GeneID == as.character(domain_control$GeneID_2[row_tmp])]
    }
  }
  domain_control = domain_control[complete.cases(domain_control),]
  #
  
  for(row_tmp in 1:dim(domain_modif)[1]){
    if(as.character(domain_modif$GeneID_1[row_tmp]) %in% specie_tmp$GeneID & as.character(domain_modif$GeneID_2[row_tmp]) %in% specie_tmp$GeneID){
      domain_modif$exp_1[row_tmp] = specie_tmp$TPM[specie_tmp$GeneID == as.character(domain_modif$GeneID_1[row_tmp])]
      domain_modif$exp_2[row_tmp] = specie_tmp$TPM[specie_tmp$GeneID == as.character(domain_modif$GeneID_2[row_tmp])]
    }
  }
  domain_modif = domain_modif[complete.cases(domain_modif),]
  domain_modif_all = domain_modif
  domain_modif_ref = domain_modif[as.character(domain_modif$GeneID_1) %in% reference_gene$GeneID,]
  
  smoothScatter(domain_control$exp_1, domain_control$exp_2,
                xlab = '',
                ylab = '',
                main = paste0(gsub('[[:digit:]]', '', domain_control$GeneID_1[1]) ,' paralogs - ',tissue_list[tissue_tmp], ' - control'), cex.main = 0.6, cex.lab = 0.6, xlim = c(0, 20), ylim = c(0, 20))
  title(xlab = paste0('reference paralog', ' expression level (log transformed)'), cex.lab = 0.4, line = 2)
  title(ylab = paste0('other paralog', ' expression level (log transformed)'), cex.lab = 0.4, line = 2)
  value_tmp1 = as.numeric(domain_control$exp_1)
  value_tmp2 = as.numeric(domain_control$exp_2)
  abline(lm(value_tmp2 ~ value_tmp1), col = 'red')
  linear_param = lm(value_tmp2 ~ value_tmp1)$coefficients
  cor_value = cor.test(value_tmp2, value_tmp1)$estimate[[1]]
  text(10, 5, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(domain_modif_ref$exp_1, domain_modif_ref$exp_2,
                xlab = '',
                ylab = '',
                main = paste0(gsub('[[:digit:]]', '', domain_control$GeneID_1[1]) ,' paralogs - ',tissue_list[tissue_tmp], ' - domain modification (only reference gene)'), cex.main = 0.6, cex.lab = 0.6, xlim = c(0, 20), ylim = c(0, 20))
  title(xlab = paste0('reference longest paralog', ' expression level (log transformed)'), cex.lab =0.4, line = 2)
  title(ylab = paste0('other paralog', ' expression level (log transformed)'), cex.lab =0.4, line = 2)
  value_tmp1 = as.numeric(domain_modif_ref$exp_1)
  value_tmp2 = as.numeric(domain_modif_ref$exp_2)
  abline(lm(value_tmp2 ~ value_tmp1), col = 'red')
  linear_param = lm(value_tmp2 ~ value_tmp1)$coefficients
  cor_value = cor.test(value_tmp2, value_tmp1)$estimate[[1]]
  text(10, 5, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(domain_modif_all$exp_1, domain_modif_all$exp_2,
                xlab = '',
                ylab = '',
                main = paste0(gsub('[[:digit:]]', '', domain_control$GeneID_1[1]) ,' paralogs - ',tissue_list[tissue_tmp], ' - domain modification (all modified gene)'), cex.main = 0.6, cex.lab = 0.6, xlim = c(0, 20), ylim = c(0, 20))
  title(xlab = paste0('longest paralog', ' expression level (log transformed)'), cex.lab =0.4, line = 2)
  title(ylab = paste0('other paralog', ' expression level (log transformed)'), cex.lab =0.4, line = 2)
  value_tmp1 = as.numeric(domain_modif_all$exp_1)
  value_tmp2 = as.numeric(domain_modif_all$exp_2)
  abline(lm(value_tmp2 ~ value_tmp1), col = 'red')
  linear_param = lm(value_tmp2 ~ value_tmp1)$coefficients
  cor_value = cor.test(value_tmp2, value_tmp1)$estimate[[1]]
  text(10, 5, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
}
dev.off()
}
}
####