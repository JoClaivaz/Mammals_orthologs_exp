'''
Joaquim Claivaz
171020

get cor.diff values for any considered dataset
'''
##FUN
#taken from rg255
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

sort_paralog = function(paralog_dataset, expression_dataset){
  
  expression_para = read.table(expression_dataset, sep = '\t', header = T)
  expression_para = aggregate(expression_para$TPM ~ expression_para$GeneID, FUN = max)
  names(expression_para) = c('GeneID', 'maxTPM')
  
  expression_para$GeneID = as.character(expression_para$GeneID)
  paralog_dataset$GeneID_1 = as.character(paralog_dataset$GeneID_1)
  paralog_dataset$GeneID_2 = as.character(paralog_dataset$GeneID_2)
  paralog_dataset$shift = as.character(paralog_dataset$shift)
  paralog_dataset$longer_domain_specie = as.character(paralog_dataset$longer_domain_specie)
  paralog_dataset$TspecF_1 = as.character(paralog_dataset$TspecF_1)
  paralog_dataset$TspecF_2 = as.character(paralog_dataset$TspecF_2)
  paralog_dataset$spec_tissue_1 = as.character(paralog_dataset$spec_tissue_1)
  paralog_dataset$spec_tissue_2 = as.character(paralog_dataset$spec_tissue_2)
  
  for (para_pair in 1:dim(paralog_dataset)[1]){
    if(expression_para$maxTPM[expression_para$GeneID == paralog_dataset$GeneID_1[para_pair]] 
       < expression_para$maxTPM[expression_para$GeneID == paralog_dataset$GeneID_2[para_pair]]
       & paralog_dataset$len_2[para_pair] >= paralog_dataset$len_1[para_pair] 
       | paralog_dataset$len_2[para_pair] > paralog_dataset$len_1[para_pair]){
      
      gene_1_tmp = paralog_dataset$GeneID_2[para_pair]
      gene_2_tmp = paralog_dataset$GeneID_1[para_pair]
      paralog_dataset$GeneID_2[para_pair] = gene_2_tmp
      paralog_dataset$GeneID_1[para_pair] = gene_1_tmp
      
      gene_1_tmp = paralog_dataset$len_2[para_pair]
      gene_2_tmp = paralog_dataset$len_1[para_pair]
      paralog_dataset$len_1[para_pair] = gene_1_tmp
      paralog_dataset$len_2[para_pair] = gene_2_tmp
      
      gene_2_tmp = paralog_dataset$tspec_1[para_pair]
      gene_1_tmp = paralog_dataset$tspec_2[para_pair]
      paralog_dataset$tspec_1[para_pair] = gene_1_tmp
      paralog_dataset$tspec_2[para_pair] = gene_2_tmp
      
      gene_2_tmp = paralog_dataset$TspecF_1[para_pair]
      gene_1_tmp = paralog_dataset$TspecF_2[para_pair]
      paralog_dataset$TspecF_1[para_pair] = gene_1_tmp
      paralog_dataset$TspecF_2[para_pair] = gene_2_tmp
      
      gene_2_tmp = paralog_dataset$spec_tissue_1[para_pair]
      gene_1_tmp = paralog_dataset$spec_tissue_2[para_pair]
      paralog_dataset$spec_tissue_1[para_pair] = gene_1_tmp
      paralog_dataset$spec_tissue_2[para_pair] = gene_2_tmp
      
      if(paralog_dataset$shift[para_pair] == 'specific_to_ubiquitous'){
        paralog_dataset$shift[para_pair] = 'ubiquitous_to_specific'
      }else if(paralog_dataset$shift[para_pair] == 'ubiquitous_to_specific'){
        paralog_dataset$shift[para_pair] = 'specific_to_ubiquitous'
        
      }
      if(paralog_dataset$longer_domain_specie[para_pair] == 'specie1'){
        paralog_dataset$longer_domain_specie[para_pair] == 'specie2'
      }else if(paralog_dataset$longer_domain_specie[para_pair] == 'specie2'){
        paralog_dataset$longer_domain_specie[para_pair] = 'specie1'
        
      }
    }
  }
  
  paralog_dataset$GeneID_1 = as.factor(paralog_dataset$GeneID_1)
  paralog_dataset$GeneID_2 = as.factor(paralog_dataset$GeneID_2)
  paralog_dataset$shift = as.factor(paralog_dataset$shift)
  paralog_dataset$longer_domain_specie = as.factor(paralog_dataset$longer_domain_specie)
  paralog_dataset$TspecF_1 = as.factor(paralog_dataset$TspecF_1)
  paralog_dataset$TspecF_2 = as.factor(paralog_dataset$TspecF_2)
  paralog_dataset$spec_tissue_1 = as.factor(paralog_dataset$spec_tissue_1)
  paralog_dataset$spec_tissue_2 = as.factor(paralog_dataset$spec_tissue_2)
  
  return(paralog_dataset)
  
}

most_occurence_vector = function(data){
  return(names(sort(table(data), decreasing = T)[1]))
}

output_cordif_test_ortho = function(specie1, specie2,
                                    regexp_considered = '_ortho_notfemale_dataset',
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
                                   all_modif = FALSE,
                                   path_folder = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/'){
  
  species_data = read.csv(paste0(path_folder, specie,
                                   regexp_considered))
  if(all_modif == FALSE){
    species_data = sort_paralog(paralog_dataset = species_data, expression_dataset = paste0('D:/UNIL/Master/Master_Project/Data/Bgee/', gsub('[[:digit:]]', '', species_data$GeneID_1[1]), '_expression_parsed'))
    paralog_reference = aggregate(species_data$GeneID_1 ~ species_data$ParalogGroup, FUN = most_occurence_vector)
    keep_gene_1 = species_data$GeneID_1 %in% paralog_reference$`species_data$GeneID_1`
    species_data = species_data[keep_gene_1,]
  }else{
    species_data = sort_paralog(paralog_dataset = species_data, expression_dataset = paste0('D:/UNIL/Master/Master_Project/Data/Bgee/', gsub('[[:digit:]]', '', species_data$GeneID_1[1]), '_expression_parsed'))
    paralog_reference = aggregate(species_data$GeneID_1 ~ species_data$ParalogGroup, FUN = most_occurence_vector)
    species_data_modif = species_data[species_data$status != 'control',]
    species_data_c = species_data[species_data$status == 'control',]
    keep_gene_1 = species_data_c$GeneID_1 %in% paralog_reference$`species_data$GeneID_1`
    species_data_c = species_data_c[keep_gene_1,]
    species_data = rbind(species_data_c, species_data_modif)
  }
  
    
  return(cor.diff.test(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
                       species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control']))
}

#

###ortholog analysis
##Run test in considered dataset
output_cordif_test_ortho(specie1 = 'HUMAN', specie2 = 'BOVIN')

###paralog analysis
##Run test in considered dataset
output_cordif_test_para(specie = 'HUMAN', 
                        regexp_considered = '_para_notfemale_dataset',
                        all_modif = FALSE)
output_cordif_test_para(specie = 'HUMAN', 
                        regexp_considered = '_para_notfemale_dataset',
                        all_modif = TRUE)
