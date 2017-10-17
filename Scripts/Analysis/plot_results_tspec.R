'''
Joaquim Claivaz
171016

Tspec analysis plot results
'''
####Ortholog plots / tspec correlation####
path_folder = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/'
different_species = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN')
central_species = c('BOVIN', 'MOUSE', 'HUMAN')
regexp_list = c('_ortho_notfemale_dataset', '_ortho_onlymale_dataset',
                '_ortho_notfemale_onlybrain_dataset', '_ortho_onlymale_onlybrain_dataset',
                '_ortho_notfemale_nottestis_dataset', '_ortho_onlymale_nottestis_dataset',
                '_ortho_notfemale_onlybrain_nottestis_dataset', '_ortho_onlymale_onlybrain_nottestis_dataset',
                '_ortho_notfemale_onlypref_nottestis_dataset', '_ortho_onlymale_onlypref_nottestis_dataset',
                '_ortho_notfemale_onlypref_dataset', '_ortho_onlymale_onlypref_dataset') 


for (regexp_out in 1:length(regexp_list)){
  pdf(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/result_cor', regexp_list[regexp_out] ,'.pdf'))
  par(mfrow = c(4,4), mai=c(0.4,0.35,0.3,0.01))
  
  for (sp1 in 1:length(central_species)){
  
    for (sp2 in 1:length(different_species)){
      
      if (central_species[sp1] != different_species[sp2]){
        options(show.error.messages = FALSE)
        try_test = try(read.csv(paste0(path_folder, central_species[sp1], '_', different_species[sp2],
                                       regexp_list[regexp_out])))
        options(show.error.messages = TRUE)
        
        if(class(try_test) != 'try-error'){
          species_data = read.csv(paste0(path_folder, central_species[sp1], '_', different_species[sp2],
                                             regexp_list[regexp_out]))
          
        
        }else{
          species_data = read.csv(paste0(path_folder, different_species[sp2], '_', central_species[sp1],
                                         regexp_list[regexp_out]))
        }
        
        species_data$X = NULL
          
        if(gsub('[[:digit:]]' ,'' , species_data$GeneID_2[1]) == central_species[sp1]){
          smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
                        xlab = '',
                        ylab = '',
                        main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
          linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
          cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
          smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                        xlab = '',
                        ylab = '',
                        main = 'Domain control group', cex.main = 0.8, cex.lab = 0.6)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
          linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
          cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
        
        }else{
          smoothScatter(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'],
                        xlab = '',
                        ylab = '',
                        main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_1[species_data$status != 'control'] ~ species_data$tspec_2[species_data$status != 'control']), col = 'red')
          linear_param = lm(species_data$tspec_1[species_data$status != 'control'] ~ species_data$tspec_2[species_data$status != 'control'])$coefficients
          cor_value = cor.test(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
          smoothScatter(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'],
                        xlab = '',
                        ylab = '',
                        main = 'Domain control group', cex.main = 0.8, cex.lab = 0.6)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_1[species_data$status == 'control'] ~ species_data$tspec_2[species_data$status == 'control']), col = 'red')
          linear_param = lm(species_data$tspec_1[species_data$status == 'control'] ~ species_data$tspec_2[species_data$status == 'control'])$coefficients
          cor_value = cor.test(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
        }
      }
    }
  }
  dev.off()
}

####Ortholog plots / correlation function of modification position####
path_folder = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/'
different_species = c('BOVIN', 'MACMU', 'MOUSE', 'PIGXX', 'HUMAN')
central_species = c('BOVIN', 'MOUSE', 'HUMAN')
regexp_list = c('_ortho_notfemale_dataset', 
                '_ortho_notfemale_onlybrain_dataset',
                '_ortho_notfemale_nottestis_dataset',
                '_ortho_notfemale_onlybrain_nottestis_dataset',
                '_ortho_notfemale_onlypref_dataset',
                '_ortho_notfemale_onlypref_nottestis_dataset')

for (regexp_out in 1:length(regexp_list)){
  pdf(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/result_modif', regexp_list[regexp_out] ,'.pdf'))
  par(mfrow = c(4,4), mai=c(0.4,0.35,0.3,0.01))
  
  for (sp1 in 1:length(central_species)){
    
    for (sp2 in 1:length(different_species)){
      
      if (central_species[sp1] != different_species[sp2]){
        options(show.error.messages = FALSE)
        try_test = try(read.csv(paste0(path_folder, central_species[sp1], '_', different_species[sp2],
                                       regexp_list[regexp_out])))
        options(show.error.messages = TRUE)
        
        if(class(try_test) != 'try-error'){
          species_data = read.csv(paste0(path_folder, central_species[sp1], '_', different_species[sp2],
                                         regexp_list[regexp_out]))
          
          
        }else{
          species_data = read.csv(paste0(path_folder, different_species[sp2], '_', central_species[sp1],
                                         regexp_list[regexp_out]))
        }
        
        species_data$X = NULL
        
        if(gsub('[[:digit:]]' ,'' , species_data$GeneID_2[1]) == central_species[sp1]){
          
          smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                        species_data$tspec_2[species_data$status == 'f-1'],
                        xlab = '',
                        ylab = '',
                        main = paste0('Domain modification group: f-1'), cex.main = 0.8)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
          linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                              species_data$tspec_1[species_data$status == 'f-1'])$coefficients
          cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                               species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
          smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                        species_data$tspec_2[species_data$status == 'int-1'],
                        xlab = '',
                        ylab = '',
                        main = paste0('Domain modification group: int-1'), cex.main = 0.8)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
          linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                              species_data$tspec_1[species_data$status == 'int-1'])$coefficients
          cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                               species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
          smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                        species_data$tspec_2[species_data$status == 'b-1'],
                        xlab = '',
                        ylab = '',
                        main = paste0('Domain modification group: b-1'), cex.main = 0.8)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
          linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                              species_data$tspec_1[species_data$status == 'b-1'])$coefficients
          cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                               species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
          smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                        xlab = '',
                        ylab = '',
                        main = 'Domain control group', cex.main = 0.8, cex.lab = 0.6)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
          linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
          cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
          
        }else{
          smoothScatter(species_data$tspec_2[species_data$status == 'f-1'],
                        species_data$tspec_1[species_data$status == 'f-1'],
                        xlab = '',
                        ylab = '',
                        main = paste0('Domain modification group: f-1'), cex.main = 0.8)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_1[species_data$status == 'f-1'] ~
                      species_data$tspec_2[species_data$status == 'f-1']), col = 'red')
          linear_param = lm(species_data$tspec_1[species_data$status == 'f-1'] ~
                              species_data$tspec_2[species_data$status == 'f-1'])$coefficients
          cor_value = cor.test(species_data$tspec_1[species_data$status == 'f-1'],
                               species_data$tspec_2[species_data$status == 'f-1'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
          smoothScatter(species_data$tspec_2[species_data$status == 'int-1'],
                        species_data$tspec_1[species_data$status == 'int-1'],
                        xlab = '',
                        ylab = '',
                        main = paste0('Domain modification group: int-1'), cex.main = 0.8)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_1[species_data$status == 'int-1'] ~
                      species_data$tspec_2[species_data$status == 'int-1']), col = 'red')
          linear_param = lm(species_data$tspec_1[species_data$status == 'int-1'] ~
                              species_data$tspec_2[species_data$status == 'int-1'])$coefficients
          cor_value = cor.test(species_data$tspec_1[species_data$status == 'int-1'],
                               species_data$tspec_2[species_data$status == 'int-1'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
          smoothScatter(species_data$tspec_2[species_data$status == 'b-1'],
                        species_data$tspec_1[species_data$status == 'b-1'],
                        xlab = '',
                        ylab = '',
                        main = paste0('Domain modification group: b-1'), cex.main = 0.8)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_1[species_data$status == 'b-1'] ~
                      species_data$tspec_2[species_data$status == 'b-1']), col = 'red')
          linear_param = lm(species_data$tspec_1[species_data$status == 'b-1'] ~
                              species_data$tspec_2[species_data$status == 'b-1'])$coefficients
          cor_value = cor.test(species_data$tspec_1[species_data$status == 'b-1'],
                               species_data$tspec_2[species_data$status == 'b-1'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
          smoothScatter(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'],
                        xlab = '',
                        ylab = '',
                        main = 'Domain control group', cex.main = 0.8, cex.lab = 0.6)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          abline(lm(species_data$tspec_1[species_data$status == 'control'] ~ species_data$tspec_2[species_data$status == 'control']), col = 'red')
          linear_param = lm(species_data$tspec_1[species_data$status == 'control'] ~ species_data$tspec_2[species_data$status == 'control'])$coefficients
          cor_value = cor.test(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'])$estimate[[1]]
          text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
          
        }
      }
    }
  }
  dev.off()
}

####Paralog / correlation Tspec####
###FUN
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
#

path_folder = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/'
central_species = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN')
regexp_list = c('_para_notfemale_dataset', '_para_onlymale_dataset',
                '_para_notfemale_onlybrain_dataset', '_para_onlymale_onlybrain_dataset',
                '_para_notfemale_nottestis_dataset', '_para_onlymale_nottestis_dataset',
                '_para_notfemale_onlybrain_nottestis_dataset', '_para_onlymale_onlybrain_nottestis_dataset',
                '_para_notfemale_onlypref_nottestis_dataset', '_para_onlymale_onlypref_nottestis_dataset',
                '_para_notfemale_onlypref_dataset', '_para_onlymale_onlypref_dataset') 

for (regexp_out in 1:length(regexp_list)){
  pdf(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/result_cor', regexp_list[regexp_out] ,'.pdf'))
  par(mfrow = c(4,4), mai=c(0.4,0.35,0.3,0.01))
  
  for (sp1 in 1:length(central_species)){
  
    species_data = read.csv(paste0(path_folder, central_species[sp1],
                                     regexp_list[regexp_out]))
    species_data$X = NULL
    
    ###Determine maximal paralog expression for a group family, according to Kryuchkova et al., 2016
    #maximal expression (reference, maximal in one state) gene is the GeneID_1
    #and longer domain is also GeneID_1 in modified gene
    #unique paralog consideration
    species_data = sort_paralog(paralog_dataset = species_data, expression_dataset = paste0('D:/UNIL/Master/Master_Project/Data/Bgee/', gsub('[[:digit:]]', '', species_data$GeneID_1[1]), '_expression_parsed'))
    paralog_reference = aggregate(species_data$GeneID_1 ~ species_data$ParalogGroup, FUN = most_occurence_vector)
    keep_gene_1 = species_data$GeneID_1 %in% paralog_reference$`species_data$GeneID_1`
    species_data = species_data[keep_gene_1,]
    #
      
    smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
                  xlab = '',
                  ylab = '',
                  main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
    title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
    title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
    abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
    linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
    cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
    text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
    
    smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                  xlab = '',
                  ylab = '',
                  main = 'Domain control group', cex.main = 0.8, cex.lab = 0.6)
    title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
    title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
    abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
    linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
    cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
    text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  }
  dev.off()
}

####Paralog plots / correlation function of modification position####
path_folder = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/'
central_species = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN')
regexp_list = c('_para_notfemale_dataset', '_para_onlymale_dataset',
                '_para_notfemale_onlybrain_dataset', '_para_onlymale_onlybrain_dataset',
                '_para_notfemale_nottestis_dataset', '_para_onlymale_nottestis_dataset',
                '_para_notfemale_onlybrain_nottestis_dataset', '_para_onlymale_onlybrain_nottestis_dataset',
                '_para_notfemale_onlypref_nottestis_dataset', '_para_onlymale_onlypref_nottestis_dataset',
                '_para_notfemale_onlypref_dataset', '_para_onlymale_onlypref_dataset') 

for (regexp_out in 1:length(regexp_list)){
  pdf(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/result_modif', regexp_list[regexp_out] ,'.pdf'))
  par(mfrow = c(4,4), mai=c(0.4,0.35,0.3,0.01))
  
  for (sp1 in 1:length(central_species)){
    
    species_data = read.csv(paste0(path_folder, central_species[sp1],
                                   regexp_list[regexp_out]))
    species_data$X = NULL
    
    ###Determine maximal paralog expression for a group family, according to Kryuchkova et al., 2016
    #maximal expression (reference, maximal in one state) gene is the GeneID_1
    #and longer domain is also GeneID_1 in modified gene
    #unique paralog consideration
    species_data = sort_paralog(paralog_dataset = species_data, expression_dataset = paste0('D:/UNIL/Master/Master_Project/Data/Bgee/', gsub('[[:digit:]]', '', species_data$GeneID_1[1]), '_expression_parsed'))
    paralog_reference = aggregate(species_data$GeneID_1 ~ species_data$ParalogGroup, FUN = most_occurence_vector)
    keep_gene_1 = species_data$GeneID_1 %in% paralog_reference$`species_data$GeneID_1`
    species_data = species_data[keep_gene_1,]
    #
    
    smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                  species_data$tspec_2[species_data$status == 'f-1'],
                  xlab = '',
                  ylab = '',
                  main = paste0('Domain modification group: f-1'), cex.main = 0.8)
    title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
    title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
    abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
    linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                        species_data$tspec_1[species_data$status == 'f-1'])$coefficients
    cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                         species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
    text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
    
    smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                  species_data$tspec_2[species_data$status == 'int-1'],
                  xlab = '',
                  ylab = '',
                  main = paste0('Domain modification group: int-1'), cex.main = 0.8)
    title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
    title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
    abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
    linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                        species_data$tspec_1[species_data$status == 'int-1'])$coefficients
    cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                         species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
    text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
    
    smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                  species_data$tspec_2[species_data$status == 'b-1'],
                  xlab = '',
                  ylab = '',
                  main = paste0('Domain modification group: b-1'), cex.main = 0.8)
    title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
    title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
    abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
    linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                        species_data$tspec_1[species_data$status == 'b-1'])$coefficients
    cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                         species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
    text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
    
    smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                  xlab = '',
                  ylab = '',
                  main = 'Domain control group', cex.main = 0.8, cex.lab = 0.6)
    title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
    title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
    abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
    linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
    cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
    text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  }
  dev.off()
}


