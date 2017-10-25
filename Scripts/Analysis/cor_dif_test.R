'''
Joaquim Claivaz
171020

get cor.diff values for any considered dataset
'''
####FUN####
#taken from rg255
{
cor.diff.test = function(x1, x2, y1, y2, method="pearson") {
  cor1 = cor.test(x1, x2, method=method)
  cor2 = cor.test(y1, y2, method=method)
  
  r1 = cor1$estimate
  r2 = cor2$estimate
  n1 = sum(complete.cases(x1, x2))
  n2 = sum(complete.cases(y1, y2))
  fisher = ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n1-3))+(1/(n2-3)))^0.5
  
  p.value = (2*(1-pnorm(abs(fisher))))
  
  result= list(
    "cor1" = list(
      "estimate" = as.numeric(cor1$estimate),
      "p.value" = cor1$p.value,
      "n" = n1
    ),
    "cor2" = list(
      "estimate" = as.numeric(cor2$estimate),
      "p.value" = cor2$p.value,
      "n" = n2
    ),
    "p.value.twosided" = as.numeric(p.value),
    "p.value.onesided" = as.numeric(p.value) / 2
  )
  cat(paste(sep="",
            "cor1: r=", format(result$cor1$estimate, digits=3), ", p=", format(result$cor1$p.value, digits=3), ", n=", result$cor1$n, "\n",
            "cor2: r=", format(result$cor2$estimate, digits=3), ", p=", format(result$cor2$p.value, digits=3), ", n=", result$cor2$n, "\n",
            "diffence: p(one-sided)=", format(result$p.value.onesided, digits=3), ", p(two-sided)=", format(result$p.value.twosided, digits=3), "\n"
  ))
  return(result);
}

#

output_cordif_test_ortho = function(specie1, specie2,
                                    regexp_considered = '_ortho_notfemale_nottestis_dataset',
                                    path_folder = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/'){
  options(show.error.messages = FALSE)
  try_test = try(read.csv(paste0(path_folder, specie1, '_', specie2,
                                 regexp_considered)))
  options(show.error.messages = TRUE)
  
  if(class(try_test) != 'try-error'){
    species_data = read.csv(paste0(path_folder, specie1, '_', specie2,
                                   regexp_considered))
    
    
  }else{
    species_data = read.csv(paste0(path_folder, specie2, '_', specie1,
                                   regexp_considered))
  }
  
  return(cor.diff.test(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
                       species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control']))
}

output_cordif_test_para = function(specie,
                                   regexp_considered,
                                   path_folder = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/'){
  
  species_data = read.csv(paste0(path_folder, specie,
                                   regexp_considered))
    
  return(cor.diff.test(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
                       species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control']))
}

number_tissue_ortholog = function(considered_species_name_1, 
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
  
  return(length(unique(as.character(specie2_data$Tissue))))
}  

number_tissue_paralog = function(considered_species_name,
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
  return(length(unique(as.character(specie_data$Tissue))))
}  
####
}
###list species available
species_list = c('BOVIN', 'GORGO', 'HUMAN', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO')

#

###ortholog analysis
##Run test in considered dataset
species_done = c()
for (sp1 in 1:length(species_list)){
  species_done = c(species_done, species_list[sp1])
  for (sp2 in 1:length(species_list)){
    if(!(species_list[sp2] %in% species_done)){
      
      options(show.error.messages = FALSE)
      try_test = try(output_cordif_test_ortho(specie1 = species_list[sp1], specie2 = species_list[sp2],
                                              regexp_considered = '_ortho_notfemale_nottestis_dataset'))
      options(show.error.messages = TRUE)
      
      if(class(try_test) != 'try-error'){
        output_cordif_test_ortho(specie1 = species_list[sp1], specie2 = species_list[sp2],
                               regexp_considered = '_ortho_notfemale_nottestis_dataset')
      }else{
        output_cordif_test_ortho(specie1 = species_list[sp2], specie2 = species_list[sp1],
                                 regexp_considered = '_ortho_notfemale_nottestis_dataset')
      }
    }
  }
}
#

###paralog analysis
##Run test in considered dataset
for (sp1 in 1:length(species_list)){
  output_cordif_test_para(specie = species_list[sp1], 
                          regexp_considered = '_para_notfemale_nottestis_dataset')
  output_cordif_test_para(specie = species_list[sp1], 
                          regexp_considered = '_para_notfemale_nottestis_dataset')
}
#

###ortholog number of tissues
species_done = c()
for (sp1 in 1:length(species_list)){
  species_done = c(species_done, species_list[sp1])
  for (sp2 in 1:length(species_list)){
    if(!(species_list[sp2] %in% species_done)){
      
      options(show.error.messages = FALSE)
      try_test = try(read.csv(paste0('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_', species_list[sp1], '_', species_list[sp2],
                                     '_domain_nomodif')))
      options(show.error.messages = TRUE)
      
      if(class(try_test) != 'try-error'){
        number_tissue_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                       expression_data_sufix = '_expression_parsed',
                       considered_species_name_1 = species_list[sp1],
                       considered_species_name_2 = species_list[sp2],
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
                       notconsidered_anat_vector = c('kidney', 'frontal cortex', 'testis'))
      }else{
        number_tissue_ortholog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                               expression_data_sufix = '_expression_parsed',
                               considered_species_name_1 = species_list[sp2],
                               considered_species_name_2 = species_list[sp1],
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
                               notconsidered_anat_vector = c('kidney', 'frontal cortex', 'testis'))  
      }
    }
  }
}
#

###paralog number of tissues
for (sp1 in 1:length(species_list)){
  number_tissue_paralog(expression_data_path_prefix = 'D:/UNIL/Master/Master_Project/Data/Bgee/',
                      expression_data_sufix = '_expression_parsed',
                      considered_species_name = species_list[sp1],
                      domain_control_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/paralog_',
                      domain_control_sufix = '_domain_nomodif',
                      domain_modif_path_prefix = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_paralog_',
                      domain_modif_sufix = '_domain_loss',
                      notconsidered_sex_vector = c('female'),
                      notconsidered_anat_vector = c('kidney', 'frontal cortex', 'testis'),
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
                                                       "19th week post-fertilization human stage (human)"))
}
#
