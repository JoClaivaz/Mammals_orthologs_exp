'''
Joaquim Claivaz
171003, last modification 171012

Tspec inference in mammals
'''
####FUN####
##Tspec analysis##
log_transformation_tpm = function(data_frame, first_column_num = 2){
  for (i in first_column_num:dim(data_frame)[2]){
    data_frame[,i] = log2(data_frame[,i] + 0.000001)
    data_frame[,i][data_frame[,i] < 1] = 0
  }
  return(data_frame)
}

#tau calculation from Nadeza work
fTau <- function(x){
  if(all(!is.na(x)))  {
    if(min(x, na.rm=TRUE) >= 0)    {
      if(max(x)!=0)      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
    } 
  } else {
    res <- NA
  } 
  return(res)
}
#

Tspec_inference = function(expression_data, first_column_num = 2){
  expression_data = log_transformation_tpm(expression_data, first_column_num)
  expression_data$tspec = NA
  for (gene_row in 1:dim(expression_data)[1]){
    expression_data$tspec[gene_row] = fTau(expression_data[gene_row, first_column_num:(dim(expression_data)[2]-1)])
  }
  
  return(expression_data)
}

specificity_inference = function(expression_data, threshold_specificity = 0.8, first_column_num){
  column_tmp = dim(expression_data)[2]-1
  
  #status inference of tspec: ubiquitous or specific (threshold 0.8)
  expression_data$TspecF = 'ubiquitous'
  expression_data$TspecF[expression_data$tspec > threshold_specificity] = 'specific'
  expression_data$TspecF = as.factor(expression_data$TspecF)
  
  #for tissue specific gene infere which tissue specificity
  expression_data$spec_tissue = 'ubiquitous'
  
  for (gene_row in 1:dim(expression_data)[1]){
    if (expression_data$TspecF[gene_row] == 'specific'){
      expression_data$spec_tissue[gene_row] = names(expression_data[gene_row, first_column_num : column_tmp])[which(expression_data[gene_row, first_column_num : column_tmp] == max(expression_data[gene_row , first_column_num : column_tmp]))] 
    }
  }
  
  return(expression_data) 
}

data_organization_Tspec_ortholog = function(considered_species_name_1, 
                                            considered_species_name_2,
                                            considered_sex_vector = F,
                                            notconsidered_sex_vector = F,
                                            considered_anat_vector = F,
                                            notconsidered_anat_vector = F,
                                            considered_devtime_vector = F,
                                            notconsidered_devtime_vector = F,
                                            expression_data_path_prefix,
                                            expression_data_sufix,
                                            domain_control_path_prefix,
                                            domain_control_sufix,
                                            domain_modif_path_prefix,
                                            domain_modif_sufix){
  
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
  specie2_data = aggregate(specie2_data$TPM ~ specie2_data$GeneID + specie2_data$AnatomicalEntityName, FUN = mean)
  names(specie2_data) = c('GeneID', 'Tissue', 'TPM')
  specie1_data = aggregate(specie1_data$TPM ~ specie1_data$GeneID + specie1_data$AnatomicalEntityName, FUN = mean)
  names(specie1_data) = c('GeneID', 'Tissue', 'TPM')
  #
  
  #group the domain inference data
  domain_modif$V6 = NULL
  domain_status = rbind(domain_control, domain_modif)
  names(domain_status) = c('GeneID_1', 'len_1', 'GeneID_2', 'len_2', 'status')
  domain_status$status = as.character(domain_status$status)
  domain_status$status[domain_status$status == 'nomodif'] = 'control'
  #
  
  #Group expression data
  #all(unique(as.character(specie2_data$Tissue)) == unique(as.character(specie1_data$Tissue)))
  names_tissue = unique(as.character(specie2_data$Tissue))
  specie1_data = aggregate(TPM ~ GeneID, data = specie1_data, paste, collapse = '_')
  specie1_data = separate(specie1_data, col = 'TPM', sep= '_', into = names_tissue)
  specie2_data = aggregate(TPM ~ GeneID, data = specie2_data, paste, collapse = '_')
  specie2_data = separate(specie2_data, col = 'TPM', sep= '_', into = names_tissue)
  #
  
  #convert tpm into numerics
  for (col_tmp in 2:dim(specie1_data)[2]){
    specie1_data[,col_tmp] = as.numeric(specie1_data[,col_tmp])
  }
  
  for (col_tmp in 2:dim(specie2_data)[2]){
    specie2_data[,col_tmp] = as.numeric(specie2_data[,col_tmp])
  }
  #
  
  #Tspec inference with log_transformation
  specie1_data = Tspec_inference(specie1_data, first_column_num = 2)
  specie2_data = Tspec_inference(specie2_data, first_column_num = 2)
  #
  
  #specificity factor inference
  specie1_data = specificity_inference(specie1_data, first_column_num = 2)
  specie2_data = specificity_inference(specie2_data, first_column_num = 2)
  #
  
  #dataframe containing Specie1 and Specie2 data
  col_keep = c('GeneID', 'tspec', 'TspecF', 'spec_tissue')
  specie1_data = specie1_data[,names(specie1_data) %in% col_keep]
  specie2_data = specie2_data[,names(specie2_data) %in% col_keep]
  species_data = rbind(specie2_data, specie1_data)
  domain_status$tspec_1 = NA
  domain_status$TspecF_1 = NA
  domain_status$spec_tissue_1 = NA
  domain_status$tspec_2 = NA
  domain_status$TspecF_2 = NA
  domain_status$spec_tissue_2 = NA
  species_data$GeneID = as.character(species_data$GeneID)
  species_data$TspecF = as.character(species_data$TspecF)
  domain_status$GeneID_1 = as.character(domain_status$GeneID_1)
  domain_status$GeneID_2 = as.character(domain_status$GeneID_2)
  
  for (gene_row in 1:dim(domain_status)[1]){
    row_sp_1 = which(species_data$GeneID == domain_status$GeneID_1[gene_row])
    row_sp_2 = which(species_data$GeneID == domain_status$GeneID_2[gene_row])
    
    if(length(row_sp_1) > 0 & length(row_sp_2) > 0){
      domain_status$tspec_1[gene_row] = species_data$tspec[row_sp_1]
      domain_status$TspecF_1[gene_row] = species_data$TspecF[row_sp_1]
      domain_status$spec_tissue_1[gene_row] = species_data$spec_tissue[row_sp_1]
      domain_status$tspec_2[gene_row] = species_data$tspec[row_sp_2]
      domain_status$TspecF_2[gene_row] = species_data$TspecF[row_sp_2]
      domain_status$spec_tissue_2[gene_row] = species_data$spec_tissue[row_sp_2]
    }
  }
  #
  
  #keep only complete case
  domain_status = domain_status[complete.cases(domain_status),]
  #
  
  #specificity shift and longer domain specie inference 
  domain_status$shift = NA
  domain_status$longer_domain_specie = NA
  
  for (pair_row in 1:dim(domain_status)[1]){
    if (domain_status$spec_tissue_1[pair_row] == domain_status$spec_tissue_2[pair_row]){
      if(domain_status$spec_tissue_1[pair_row] != 'ubiquitous'){
        domain_status$shift[pair_row] = 'specific_no_shift'
      }else{
        domain_status$shift[pair_row] = 'ubiquitous_no_shift'
      }
    }
    else if (all(domain_status$spec_tissue_1[pair_row] != 'ubiquitous',
                 domain_status$spec_tissue_2[pair_row] != 'ubiquitous')){
      domain_status$shift[pair_row] = 'specificity_shift'
    }
    else if (domain_status$spec_tissue_1[pair_row] != 'ubiquitous'){
      domain_status$shift[pair_row] = 'specific_to_ubiquitous'
    }
    else if (domain_status$spec_tissue_2[pair_row] != 'ubiquitous'){
      domain_status$shift[pair_row] = 'ubiquitous_to_specific'
    }
    else{
      domain_status$shift[pair_row] = 'undefined'
    }
    if (domain_status$len_1[pair_row] != domain_status$len_2[pair_row]){
      if (domain_status$len_1[pair_row] > domain_status$len_2[pair_row]){
        domain_status$longer_domain_specie[pair_row] = 'specie1'
      }else{
        domain_status$longer_domain_specie[pair_row] = 'specie2'
      }
    }else{
      domain_status$longer_domain_specie[pair_row] = 'control'
    }
  }
  #
  
  domain_status$GeneID_1 = as.factor(domain_status$GeneID_1)
  domain_status$GeneID_2 = as.factor(domain_status$GeneID_2)
  domain_status$status = as.factor(domain_status$status)
  domain_status$TspecF_1 = as.factor(domain_status$TspecF_1)
  domain_status$TspecF_2 = as.factor(domain_status$TspecF_2)
  domain_status$spec_tissue_1 = as.factor(domain_status$spec_tissue_1)
  domain_status$spec_tissue_2 = as.factor(domain_status$spec_tissue_2)
  domain_status$longer_domain_specie = as.factor(domain_status$longer_domain_specie)
  domain_status$shift = as.factor(domain_status$shift)
  
  return(domain_status)
}

data_organization_Tspec_paralog = function(considered_species_name,
                                           considered_sex_vector = F,
                                           notconsidered_sex_vector = F,
                                           considered_anat_vector = F,
                                           notconsidered_anat_vector = F,
                                           considered_devtime_vector = F,
                                           notconsidered_devtime_vector = F,
                                           expression_data_path_prefix,
                                           expression_data_sufix,
                                           domain_control_path_prefix,
                                           domain_control_sufix,
                                           domain_modif_path_prefix,
                                           domain_modif_sufix){
  
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
  
  #group the domain inference data
  domain_control$V8 = NULL
  domain_status = rbind(domain_control, domain_modif)
  domain_status$V5 = as.character(domain_status$V5)
  domain_status$V5[domain_status$V5 == 'paralog'] = 'control'
  names(domain_status) = c('GeneID_1', 'len_1', 'GeneID_2', 'len_2', 'status', 'ParalogGroup', 'SpeciesParalogInference')
  #
  
  #Group expression data
  #all(unique(as.character(specie2_data$Tissue)) == unique(as.character(specie1_data$Tissue)))
  names_tissue = unique(as.character(specie_data$Tissue))
  specie_data = aggregate(TPM ~ GeneID, data = specie_data, paste, collapse = '_')
  specie_data = separate(specie_data, col = 'TPM', sep= '_', into = names_tissue)
  #
  
  #convert tpm into numerics
  for (col_tmp in 2:dim(specie_data)[2]){
    specie_data[,col_tmp] = as.numeric(specie_data[,col_tmp])
  }
  #
  
  #Tspec inference with log_transformation
  specie_data = Tspec_inference(specie_data, first_column_num = 2)
  #
  
  #specificity factor inference
  specie_data = specificity_inference(specie_data, first_column_num = 2)
  #
  
  #dataframe containing Specie1 and Specie2 data
  col_keep = c('GeneID', 'tspec', 'TspecF', 'spec_tissue')
  specie_data = specie_data[,names(specie_data) %in% col_keep]
  domain_status$tspec_1 = NA
  domain_status$TspecF_1 = NA
  domain_status$spec_tissue_1 = NA
  domain_status$tspec_2 = NA
  domain_status$TspecF_2 = NA
  domain_status$spec_tissue_2 = NA
  specie_data$GeneID = as.character(specie_data$GeneID)
  specie_data$TspecF = as.character(specie_data$TspecF)
  domain_status$GeneID_1 = as.character(domain_status$GeneID_1)
  domain_status$GeneID_2 = as.character(domain_status$GeneID_2)
  
  for (gene_row in 1:dim(domain_status)[1]){
    row_sp_1 = which(specie_data$GeneID == domain_status$GeneID_1[gene_row])
    row_sp_2 = which(specie_data$GeneID == domain_status$GeneID_2[gene_row])
    
    if(length(row_sp_1) > 0 & length(row_sp_2) > 0){
      domain_status$tspec_1[gene_row] = specie_data$tspec[row_sp_1]
      domain_status$TspecF_1[gene_row] = specie_data$TspecF[row_sp_1]
      domain_status$spec_tissue_1[gene_row] = specie_data$spec_tissue[row_sp_1]
      domain_status$tspec_2[gene_row] = specie_data$tspec[row_sp_2]
      domain_status$TspecF_2[gene_row] = specie_data$TspecF[row_sp_2]
      domain_status$spec_tissue_2[gene_row] = specie_data$spec_tissue[row_sp_2]
    }
  }
  #
  
  #keep only complete case
  domain_status = domain_status[complete.cases(domain_status),]
  #
  
  #specificity shift and longer domain specie inference 
  domain_status$shift = NA
  domain_status$longer_domain_specie = NA
  
  for (pair_row in 1:dim(domain_status)[1]){
    if (domain_status$spec_tissue_1[pair_row] == domain_status$spec_tissue_2[pair_row]){
      if(domain_status$spec_tissue_1[pair_row] != 'ubiquitous'){
        domain_status$shift[pair_row] = 'specific_no_shift'
      }else{
        domain_status$shift[pair_row] = 'ubiquitous_no_shift'
      }
    }
    else if (all(domain_status$spec_tissue_1[pair_row] != 'ubiquitous',
                 domain_status$spec_tissue_2[pair_row] != 'ubiquitous')){
      domain_status$shift[pair_row] = 'specificity_shift'
    }
    else if (domain_status$spec_tissue_1[pair_row] != 'ubiquitous'){
      domain_status$shift[pair_row] = 'specific_to_ubiquitous'
    }
    else if (domain_status$spec_tissue_2[pair_row] != 'ubiquitous'){
      domain_status$shift[pair_row] = 'ubiquitous_to_specific'
    }
    else{
      domain_status$shift[pair_row] = 'undefined'
    }
    if (domain_status$len_1[pair_row] != domain_status$len_2[pair_row]){
      if (domain_status$len_1[pair_row] > domain_status$len_2[pair_row]){
        domain_status$longer_domain_specie[pair_row] = 'specie1'
      }else{
        domain_status$longer_domain_specie[pair_row] = 'specie2'
      }
    }else{
      domain_status$longer_domain_specie[pair_row] = 'control'
    }
  }
  #
  
  domain_status$GeneID_1 = as.factor(domain_status$GeneID_1)
  domain_status$GeneID_2 = as.factor(domain_status$GeneID_2)
  domain_status$status = as.factor(domain_status$status)
  domain_status$TspecF_1 = as.factor(domain_status$TspecF_1)
  domain_status$TspecF_2 = as.factor(domain_status$TspecF_2)
  domain_status$spec_tissue_1 = as.factor(domain_status$spec_tissue_1)
  domain_status$spec_tissue_2 = as.factor(domain_status$spec_tissue_2)
  domain_status$longer_domain_specie = as.factor(domain_status$longer_domain_specie)
  domain_status$shift = as.factor(domain_status$shift)
  domain_status$ParalogGroup = as.factor(domain_status$ParalogGroup)
  
  return(domain_status)
}
##
####

####Library needed####
require(tidyr)
####

####Data organization####
# human_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/HUMAN_expression_parsed', header = T, sep = '\t')
# mouse_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/MOUSE_expression_parsed', header = T, sep = '\t')
# unique(mouse_test$AnatomicalEntityName)[unique(mouse_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(mouse_test$StageName)[unique(mouse_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(mouse_test$Sex)[unique(mouse_test$Sex) %in% unique(human_test$Sex)]

species_vector = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN') 
specie_done = c()
#not_female
for (cons_specie1 in 1:length(species_vector)){
  specie_done = c(specie_done, species_vector[cons_specie1])
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            notconsidered_sex_vector = c('female'),
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
                                            notconsidered_anat_vector = c('kidney', 'frontal cortex')),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_notfemale_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            considered_sex_vector = c('male'),
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
                                            notconsidered_anat_vector = c('kidney', 'frontal cortex')),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_onlymale_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            notconsidered_sex_vector = c('female'),
                                            notconsidered_anat_vector = c('testis', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_notfemale_nottestis_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            considered_sex_vector = c('male'),
                                            notconsidered_anat_vector = c('testis', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_onlymale_nottestis_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            notconsidered_sex_vector = c('female'),
                                            notconsidered_anat_vector = c('prefrontal cortex', 'cerebellum', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_notfemale_onlybrain_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            considered_sex_vector = c('male'),
                                            notconsidered_anat_vector = c('prefrontal cortex', 'cerebellum', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_onlymale_onlybrain_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            notconsidered_sex_vector = c('female'),
                                            notconsidered_anat_vector = c('prefrontal cortex', 'cerebellum', 'testis', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_notfemale_onlybrain_nottestis_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            considered_sex_vector = c('male'),
                                            notconsidered_anat_vector = c('prefrontal cortex', 'cerebellum', 'testis', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_onlymale_onlybrain_nottestis_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            notconsidered_sex_vector = c('female'),
                                            notconsidered_anat_vector = c('brain', 'cerebellum', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_notfemale_onlypref_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            considered_sex_vector = c('male'),
                                            notconsidered_anat_vector = c('brain', 'cerebellum', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_onlymale_onlypref_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            notconsidered_sex_vector = c('female'),
                                            notconsidered_anat_vector = c('brain', 'cerebellum', 'testis', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_notfemale_onlypref_nottestis_dataset'))
  
  write.csv(data_organization_Tspec_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                            expression_data_sufix = '_expression_parsed',
                                            considered_species_name = species_vector[cons_specie1],
                                            domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                                            domain_control_sufix = '_domain_nomodif',
                                            domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                                            domain_modif_sufix = '_domain_loss',
                                            considered_sex_vector = c('male'),
                                            notconsidered_anat_vector = c('brain', 'cerebellum', 'testis', 'kidney', 'frontal cortex'),
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
                                                                             "19th week post-fertilization human stage (human)")),
            file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                          species_vector[cons_specie1], '_', 'para_onlymale_onlypref_nottestis_dataset'))
  
  for (cons_specie2 in 1:length(species_vector)){
    if (!(species_vector[cons_specie2] %in% specie_done)){
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 notconsidered_sex_vector = c('female'),
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
                                                 notconsidered_anat_vector = c('kidney', 'frontal cortex')), 
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_notfemale_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 considered_sex_vector = c('male'),
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
                                                 notconsidered_anat_vector = c('kidney', 'frontal cortex')),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_onlymale_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 notconsidered_sex_vector = c('female'),
                                                 notconsidered_anat_vector = c('testis', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_notfemale_nottestis_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 considered_sex_vector = c('male'),
                                                 notconsidered_anat_vector = c('testis', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_onlymale_nottestis_dataset'))
      
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 notconsidered_sex_vector = c('female'),
                                                 notconsidered_anat_vector = c('prefrontal cortex', 'cerebellum', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_notfemale_onlybrain_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 considered_sex_vector = c('male'),
                                                 notconsidered_anat_vector = c('prefrontal cortex', 'cerebellum', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_onlymale_onlybrain_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 notconsidered_sex_vector = c('female'),
                                                 notconsidered_anat_vector = c('prefrontal cortex', 'cerebellum', 'testis', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_notfemale_onlybrain_nottestis_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 considered_sex_vector = c('male'),
                                                 notconsidered_anat_vector = c('prefrontal cortex', 'cerebellum', 'testis', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_onlymale_onlybrain_nottestis_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 notconsidered_sex_vector = c('female'),
                                                 notconsidered_anat_vector = c('brain', 'cerebellum', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_notfemale_onlypref_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 considered_sex_vector = c('male'),
                                                 notconsidered_anat_vector = c('brain', 'cerebellum', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_onlymale_onlypref_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 notconsidered_sex_vector = c('female'),
                                                 notconsidered_anat_vector = c('brain', 'cerebellum', 'testis', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_notfemale_onlypref_nottestis_dataset'))
      
      write.csv(data_organization_Tspec_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                                                 expression_data_sufix = '_expression_parsed',
                                                 considered_species_name_1 = species_vector[cons_specie1],
                                                 considered_species_name_2 = species_vector[cons_specie2],
                                                 domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_',
                                                 domain_control_sufix = '_domain_nomodif',
                                                 domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_',
                                                 domain_modif_sufix = '_domain_loss',
                                                 considered_sex_vector = c('male'),
                                                 notconsidered_anat_vector = c('brain', 'cerebellum', 'testis', 'kidney', 'frontal cortex'),
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
                                                                                  "19th week post-fertilization human stage (human)")),
                file = paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', 
                              species_vector[cons_specie1], '_', species_vector[cons_specie2], '_', 'ortho_onlymale_onlypref_nottestis_dataset'))
      
    }  
  }
}
####