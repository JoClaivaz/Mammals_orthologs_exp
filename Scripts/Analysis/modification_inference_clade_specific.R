'''
Joaquim Claivaz
171112

infer domain modification gene specific to a clade
only orthologs unduplicated in whole species
only gene which are present in one domain modification or control group are taken into account
'''


####Library needed####
require(tidyr)

#####Data organization#####
####HUMAN####
#store gene status (control or modification) into list for each species
list_human = list()
other_species = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO')

for (species_tmp in 1: length(other_species)){
  
  considered_species = other_species[species_tmp]
  
  data_tmp = read.table(paste0('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_HUMAN_', considered_species, '_domain_loss'), sep = '\t')
  data_tmp$status = 'modif'
  data_control = read.table(paste0('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_HUMAN_', considered_species, '_domain_nomodif'), sep = '\t')
  data_control$status = 'control'
  keep_control = c(1, 3, 6)
  keep_tmp = c(1, 3, 7)
  data_tmp = data_tmp[,keep_tmp]
  data_control = data_control[,keep_control]
  data_tmp = rbind(data_tmp, data_control)
  
  list_human[[considered_species]] = data_tmp
} 

#create dataframe from the list species
data_human = data.frame(GeneID = c(as.character(list_human[["BOVIN"]]$V1),
                                   as.character(list_human[["GORGO"]]$V3),
                                   as.character(list_human[["MACMU"]]$V3),
                                   as.character(list_human[["MONDO"]]$V1),
                                   as.character(list_human[["MOUSE"]]$V3),
                                   as.character(list_human[["PANTR"]]$V1),
                                   as.character(list_human[["PIGXX"]]$V1),
                                   as.character(list_human[["RATNO"]]$V3)),
                        status = c(as.character(list_human[['BOVIN']]$status),
                                   as.character(list_human[['GORGO']]$status),
                                   as.character(list_human[['MACMU']]$status),
                                   as.character(list_human[['MONDO']]$status),
                                   as.character(list_human[['MOUSE']]$status),
                                   as.character(list_human[['PANTR']]$status),
                                   as.character(list_human[['PIGXX']]$status),
                                   as.character(list_human[['RATNO']]$status)),
                        species = c(rep('BOVIN', length(list_human[["BOVIN"]]$V1)),
                                    rep('GORGO',length(list_human[["GORGO"]]$V3)),
                                    rep('MACMU',length(list_human[["MACMU"]]$V3)),
                                    rep('MONDO', length(list_human[["MONDO"]]$V1)),
                                    rep('MOUSE', length(list_human[["MOUSE"]]$V3)),
                                    rep('PANTR', length(list_human[["PANTR"]]$V1)),
                                    rep('PIGXX', length(list_human[["PIGXX"]]$V1)),
                                    rep('RATNO', length(list_human[["RATNO"]]$V3))))

#organize dataframe to have column with the species name and containing status
data_human = spread(data_human, key = 'species', value = 'status')
#keep only complete case which mean only gene present in 1 domain modification or control group in any species are considered
data_human = data_human[complete.cases(data_human),]
####
####MOUSE####
#store gene status (control or modification) into list for each species
list_mouse = list()
other_species = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'HUMAN', 'PANTR', 'PIGXX', 'RATNO')

for (species_tmp in 1: length(other_species)){
  
  considered_species = other_species[species_tmp]
  
  options(show.error.messages = FALSE)
  try_test = try(read.table(paste0('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_MOUSE_', considered_species, '_domain_loss'), sep = '\t'))
  options(show.error.messages = TRUE)
  
  if(class(try_test) != 'try-error'){
    data_tmp = read.table(paste0('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_MOUSE_', considered_species, '_domain_loss'), sep = '\t')
    data_control = read.table(paste0('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_MOUSE_', considered_species, '_domain_nomodif'), sep = '\t')
  }else{
    data_tmp = read.table(paste0('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/putative_ortholog_', considered_species, '_MOUSE_domain_loss'), sep = '\t')
    data_control = read.table(paste0('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_', considered_species, '_MOUSE_domain_nomodif'), sep = '\t')
  }
  data_tmp$status = 'modif'
  data_control$status = 'control'
  keep_control = c(1, 3, 6)
  keep_tmp = c(1, 3, 7)
  data_tmp = data_tmp[,keep_tmp]
  data_control = data_control[,keep_control]
  data_tmp = rbind(data_tmp, data_control)
  
  list_mouse[[considered_species]] = data_tmp
} 

#create dataframe from the list species
data_human = data.frame(GeneID = c(as.character(list_human[["BOVIN"]]$V1),
                                   as.character(list_human[["GORGO"]]$V3),
                                   as.character(list_human[["MACMU"]]$V3),
                                   as.character(list_human[["MONDO"]]$V1),
                                   as.character(list_human[["MOUSE"]]$V3),
                                   as.character(list_human[["PANTR"]]$V1),
                                   as.character(list_human[["PIGXX"]]$V1),
                                   as.character(list_human[["RATNO"]]$V3)),
                        status = c(as.character(list_human[['BOVIN']]$status),
                                   as.character(list_human[['GORGO']]$status),
                                   as.character(list_human[['MACMU']]$status),
                                   as.character(list_human[['MONDO']]$status),
                                   as.character(list_human[['MOUSE']]$status),
                                   as.character(list_human[['PANTR']]$status),
                                   as.character(list_human[['PIGXX']]$status),
                                   as.character(list_human[['RATNO']]$status)),
                        species = c(rep('BOVIN', length(list_human[["BOVIN"]]$V1)),
                                    rep('GORGO',length(list_human[["GORGO"]]$V3)),
                                    rep('MACMU',length(list_human[["MACMU"]]$V3)),
                                    rep('MONDO', length(list_human[["MONDO"]]$V1)),
                                    rep('MOUSE', length(list_human[["MOUSE"]]$V3)),
                                    rep('PANTR', length(list_human[["PANTR"]]$V1)),
                                    rep('PIGXX', length(list_human[["PIGXX"]]$V1)),
                                    rep('RATNO', length(list_human[["RATNO"]]$V3))))

data_mouse = data.frame(GeneID = c(as.character(list_mouse[["BOVIN"]]$V1),
                                   as.character(list_mouse[["GORGO"]]$V1),
                                   as.character(list_mouse[["MACMU"]]$V1),
                                   as.character(list_mouse[["MONDO"]]$V1),
                                   as.character(list_mouse[["HUMAN"]]$V1),
                                   as.character(list_mouse[["PANTR"]]$V1),
                                   as.character(list_mouse[["PIGXX"]]$V1),
                                   as.character(list_mouse[["RATNO"]]$V1)),
                        status = c(as.character(list_mouse[['BOVIN']]$status),
                                   as.character(list_mouse[['GORGO']]$status),
                                   as.character(list_mouse[['MACMU']]$status),
                                   as.character(list_mouse[['MONDO']]$status),
                                   as.character(list_mouse[['HUMAN']]$status),
                                   as.character(list_mouse[['PANTR']]$status),
                                   as.character(list_mouse[['PIGXX']]$status),
                                   as.character(list_mouse[['RATNO']]$status)),
                        species = c(rep('BOVIN', length(list_mouse[["BOVIN"]]$V1)),
                                    rep('GORGO',length(list_mouse[["GORGO"]]$V3)),
                                    rep('MACMU',length(list_mouse[["MACMU"]]$V3)),
                                    rep('MONDO', length(list_mouse[["MONDO"]]$V1)),
                                    rep('HUMAN', length(list_mouse[["HUMAN"]]$V3)),
                                    rep('PANTR', length(list_mouse[["PANTR"]]$V1)),
                                    rep('PIGXX', length(list_mouse[["PIGXX"]]$V1)),
                                    rep('RATNO', length(list_mouse[["RATNO"]]$V3))))

#organize dataframe to have column with the species name and containing status
data_human = spread(data_human, key = 'species', value = 'status')
data_mouse = spread(data_mouse, key = 'species', value = 'status')
#keep only complete case which mean only gene present in 1 domain modification or control group in any species are considered
data_human = data_human[complete.cases(data_human),]
data_mouse = data_mouse[complete.cases(data_mouse),]
####

#####Generate gene list specific to different clade#####
data_human$GeneID = as.character(data_human$GeneID)
data_mouse$GeneID = as.character(data_mouse$GeneID)

primate_gene = data_human$GeneID[data_human$BOVIN == 'modif'
                                 & data_human$MONDO == 'modif'
                                 & data_human$PIGXX == 'modif'
                                 & data_human$MOUSE == 'modif'
                                 & data_human$RATNO == 'modif'
                                 & data_human$MACMU == 'control'
                                 & data_human$GORGO == 'control'
                                 & data_human$PANTR == 'control']

hominidae_gene = data_human$GeneID[data_human$BOVIN == 'modif'
                                   & data_human$MONDO == 'modif'
                                   & data_human$PIGXX == 'modif'
                                   & data_human$MOUSE == 'modif'
                                   & data_human$RATNO == 'modif'
                                   & data_human$MACMU == 'modif'
                                   & data_human$GORGO == 'control'
                                   & data_human$PANTR == 'control']

human_gene = data_human$GeneID[data_human$BOVIN == 'modif'
                               & data_human$MONDO == 'modif'
                               & data_human$PIGXX == 'modif'
                               & data_human$MOUSE == 'modif'
                               & data_human$RATNO == 'modif'
                               & data_human$MACMU == 'modif'
                               & data_human$GORGO == 'modif'
                               & data_human$PANTR == 'modif']

muridae_gene = data_mouse$GeneID[data_mouse$BOVIN == 'modif'
                                 & data_mouse$MONDO == 'modif'
                                 & data_mouse$PIGXX == 'modif'
                                 & data_mouse$HUMAN == 'modif'
                                 & data_mouse$RATNO == 'control'
                                 & data_mouse$MACMU == 'modif'
                                 & data_mouse$GORGO == 'modif'
                                 & data_mouse$PANTR == 'modif']

mouse_gene = data_mouse$GeneID[data_mouse$BOVIN == 'modif'
                               & data_mouse$MONDO == 'modif'
                               & data_mouse$PIGXX == 'modif'
                               & data_mouse$HUMAN == 'modif'
                               & data_mouse$RATNO == 'modif'
                               & data_mouse$MACMU == 'modif'
                               & data_mouse$GORGO == 'modif'
                               & data_mouse$PANTR == 'modif']

#####Compare tissue specificity index conservation between group, human vs mouse dataset#####
####with testis####
expression_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/MOUSE_HUMAN_ortho_notfemale_dataset')
expression_data$X = NULL
keep_expresssion = c(1 ,3 , 6, 9)
expression_data = expression_data[,keep_expresssion]

expression_primate = expression_data[expression_data$GeneID_2 %in% primate_gene,]
expression_hominidae = expression_data[expression_data$GeneID_2 %in% hominidae_gene,]
expression_human = expression_data[expression_data$GeneID_2 %in% human_gene,]
expression_muridae = expression_data[expression_data$GeneID_1 %in% muridae_gene,]
expression_mouse = expression_data[expression_data$GeneID_1 %in% mouse_gene,]

expression_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/PANTR_HUMAN_ortho_notfemale_dataset')
expression_data$X = NULL
keep_expresssion = c(1 ,3 , 6, 9)
expression_data = expression_data[,keep_expresssion]
expression_human_pantr = expression_data[expression_data$GeneID_1 %in% human_gene,]

expression_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/MACMU_HUMAN_ortho_notfemale_dataset')
expression_data$X = NULL
keep_expresssion = c(1 ,3 , 6, 9)
expression_data = expression_data[,keep_expresssion]
expression_hominidae_macmu = expression_data[expression_data$GeneID_2 %in% hominidae_gene,]

expression_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/MOUSE_RATNO_ortho_notfemale_dataset')
expression_data$X = NULL
keep_expresssion = c(1 ,3 , 6, 9)
expression_data = expression_data[,keep_expresssion]
expression_mouse_ratno = expression_data[expression_data$GeneID_1 %in% mouse_gene,]

####without testis####
expression_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/MOUSE_HUMAN_ortho_notfemale_nottestis_dataset')
expression_data$X = NULL
keep_expresssion = c(1 ,3 , 6, 9)
expression_data = expression_data[,keep_expresssion]

expression_wt_primate = expression_data[expression_data$GeneID_2 %in% primate_gene,]
expression_wt_hominidae = expression_data[expression_data$GeneID_2 %in% hominidae_gene,]
expression_wt_human = expression_data[expression_data$GeneID_2 %in% human_gene,]
expression_wt_muridae = expression_data[expression_data$GeneID_1 %in% muridae_gene,]
expression_wt_mouse = expression_data[expression_data$GeneID_1 %in% mouse_gene,]

expression_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/PANTR_HUMAN_ortho_notfemale_nottestis_dataset')
expression_data$X = NULL
keep_expresssion = c(1 ,3 , 6, 9)
expression_data = expression_data[,keep_expresssion]
expression_wt_human_pantr = expression_data[expression_data$GeneID_1 %in% human_gene,]

expression_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/MACMU_HUMAN_ortho_notfemale_nottestis_dataset')
expression_data$X = NULL
keep_expresssion = c(1 ,3 , 6, 9)
expression_data = expression_data[,keep_expresssion]
expression_wt_hominidae_macmu = expression_data[expression_data$GeneID_2 %in% hominidae_gene,]

expression_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/MOUSE_RATNO_ortho_notfemale_nottestis_dataset')
expression_data$X = NULL
keep_expresssion = c(1 ,3 , 6, 9)
expression_data = expression_data[,keep_expresssion]
expression_wt_mouse_ratno = expression_data[expression_data$GeneID_1 %in% mouse_gene,]

#plot
{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'with_testis' ,'.pdf'))
  par(mfrow = c(5,4), mai=c(0.4,0.35,0.3,0.01))
  
  smoothScatter(expression_primate$tspec_1, 
                expression_primate$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Primate clade (with testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_primate$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_primate$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_primate$tspec_2 ~ 
              expression_primate$tspec_1), col = 'red')
  linear_param = lm(expression_primate$tspec_2 ~ 
                      expression_primate$tspec_1)$coefficients
  cor_value = cor.test(expression_primate$tspec_2, 
                       expression_primate$tspec_1)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_hominidae$tspec_1, 
                expression_hominidae$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Hominidae clade (with testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_hominidae$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_hominidae$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_hominidae$tspec_2 ~ 
              expression_hominidae$tspec_1), col = 'red')
  linear_param = lm(expression_hominidae$tspec_2 ~ 
                      expression_hominidae$tspec_1)$coefficients
  cor_value = cor.test(expression_hominidae$tspec_2, 
                       expression_hominidae$tspec_1)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_hominidae_macmu$tspec_1, 
                expression_hominidae_macmu$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Hominidae clade (with testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_hominidae_macmu$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_hominidae_macmu$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_hominidae_macmu$tspec_2 ~ 
              expression_hominidae_macmu$tspec_1), col = 'red')
  linear_param = lm(expression_hominidae_macmu$tspec_2 ~ 
                      expression_hominidae_macmu$tspec_1)$coefficients
  cor_value = cor.test(expression_hominidae_macmu$tspec_2, 
                       expression_hominidae_macmu$tspec_1)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_human$tspec_1, 
                expression_human$tspec_2,
                xlab = '',
                ylab = '',
                main = 'HUMAN species (with testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_human$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_human$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_human$tspec_2 ~ 
              expression_human$tspec_1), col = 'red')
  linear_param = lm(expression_human$tspec_2 ~ 
                      expression_human$tspec_1)$coefficients
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n'), cex = 0.8, col = 'red')
  
  smoothScatter(expression_human_pantr$tspec_2, 
                expression_human_pantr$tspec_1,
                xlab = '',
                ylab = '',
                main = 'HUMAN species (with testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_human_pantr$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_human_pantr$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_human_pantr$tspec_1 ~ 
              expression_human_pantr$tspec_2), col = 'red')
  linear_param = lm(expression_human_pantr$tspec_1 ~ 
                      expression_human_pantr$tspec_2)$coefficients
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n'), cex = 0.8, col = 'red')
  
  smoothScatter(expression_muridae$tspec_2, 
                expression_muridae$tspec_1,
                xlab = '',
                ylab = '',
                main = 'Muridae clade (with testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_muridae$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_muridae$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_muridae$tspec_1 ~ 
              expression_muridae$tspec_2), col = 'red')
  linear_param = lm(expression_muridae$tspec_1 ~ 
                      expression_muridae$tspec_2)$coefficients
  cor_value = cor.test(expression_muridae$tspec_1, 
                       expression_muridae$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_mouse$tspec_2, 
                expression_mouse$tspec_1,
                xlab = '',
                ylab = '',
                main = 'MOUSE species (with testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_mouse$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_mouse$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_mouse$tspec_1 ~ 
              expression_mouse$tspec_2), col = 'red')
  linear_param = lm(expression_mouse$tspec_1 ~ 
                      expression_mouse$tspec_2)$coefficients
  cor_value = cor.test(expression_mouse$tspec_1, 
                       expression_mouse$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_mouse_ratno$tspec_2, 
                expression_mouse_ratno$tspec_1,
                xlab = '',
                ylab = '',
                main = 'MOUSE species (with testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_mouse_ratno$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_mouse_ratno$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_mouse_ratno$tspec_1 ~ 
              expression_mouse_ratno$tspec_2), col = 'red')
  linear_param = lm(expression_mouse_ratno$tspec_1 ~ 
                      expression_mouse_ratno$tspec_2)$coefficients
  cor_value = cor.test(expression_mouse_ratno$tspec_1, 
                       expression_mouse_ratno$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_wt_primate$tspec_1, 
                expression_wt_primate$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Primate clade (without testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_wt_primate$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_wt_primate$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_wt_primate$tspec_2 ~ 
              expression_wt_primate$tspec_1), col = 'red')
  linear_param = lm(expression_wt_primate$tspec_2 ~ 
                      expression_wt_primate$tspec_1)$coefficients
  cor_value = cor.test(expression_wt_primate$tspec_2, 
                       expression_wt_primate$tspec_1)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_wt_hominidae$tspec_1, 
                expression_wt_hominidae$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Hominidae clade (without testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_wt_hominidae$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_wt_hominidae$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_wt_hominidae$tspec_2 ~ 
              expression_wt_hominidae$tspec_1), col = 'red')
  linear_param = lm(expression_wt_hominidae$tspec_2 ~ 
                      expression_wt_hominidae$tspec_1)$coefficients
  cor_value = cor.test(expression_wt_hominidae$tspec_2, 
                       expression_wt_hominidae$tspec_1)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_wt_hominidae_macmu$tspec_1, 
                expression_wt_hominidae_macmu$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Hominidae clade (without testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_wt_hominidae_macmu$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_wt_hominidae_macmu$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_wt_hominidae_macmu$tspec_2 ~ 
              expression_wt_hominidae_macmu$tspec_1), col = 'red')
  linear_param = lm(expression_wt_hominidae_macmu$tspec_2 ~ 
                      expression_wt_hominidae_macmu$tspec_1)$coefficients
  cor_value = cor.test(expression_wt_hominidae_macmu$tspec_2, 
                       expression_wt_hominidae_macmu$tspec_1)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_wt_human$tspec_1, 
                expression_wt_human$tspec_2,
                xlab = '',
                ylab = '',
                main = 'HUMAN species (without testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_wt_human$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_wt_human$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_wt_human$tspec_2 ~ 
              expression_wt_human$tspec_1), col = 'red')
  linear_param = lm(expression_wt_human$tspec_2 ~ 
                      expression_wt_human$tspec_1)$coefficients
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n'), cex = 0.8, col = 'red')
  
  smoothScatter(expression_wt_human_pantr$tspec_2, 
                expression_wt_human_pantr$tspec_1,
                xlab = '',
                ylab = '',
                main = 'HUMAN species (without testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_wt_human_pantr$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_wt_human_pantr$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_wt_human_pantr$tspec_1 ~ 
              expression_wt_human_pantr$tspec_2), col = 'red')
  linear_param = lm(expression_wt_human_pantr$tspec_1 ~ 
                      expression_wt_human_pantr$tspec_2)$coefficients
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n'), cex = 0.8, col = 'red')
  
  smoothScatter(expression_wt_muridae$tspec_2, 
                expression_wt_muridae$tspec_1,
                xlab = '',
                ylab = '',
                main = 'Muridae clade (without testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_wt_muridae$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_wt_muridae$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_wt_muridae$tspec_1 ~ 
              expression_wt_muridae$tspec_2), col = 'red')
  linear_param = lm(expression_wt_muridae$tspec_1 ~ 
                      expression_wt_muridae$tspec_2)$coefficients
  cor_value = cor.test(expression_wt_muridae$tspec_1, 
                       expression_wt_muridae$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_wt_mouse$tspec_2, 
                expression_wt_mouse$tspec_1,
                xlab = '',
                ylab = '',
                main = 'MOUSE species (without testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_wt_mouse$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_wt_mouse$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_wt_mouse$tspec_1 ~ 
              expression_wt_mouse$tspec_2), col = 'red')
  linear_param = lm(expression_wt_mouse$tspec_1 ~ 
                      expression_wt_mouse$tspec_2)$coefficients
  cor_value = cor.test(expression_wt_mouse$tspec_1, 
                       expression_wt_mouse$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(expression_wt_mouse_ratno$tspec_2, 
                expression_wt_mouse_ratno$tspec_1,
                xlab = '',
                ylab = '',
                main = 'MOUSE species (without testis)', cex.main = 0.8, cex.lab = 0.6, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,expression_wt_mouse_ratno$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,expression_wt_mouse_ratno$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(expression_wt_mouse_ratno$tspec_1 ~ 
              expression_wt_mouse_ratno$tspec_2), col = 'red')
  linear_param = lm(expression_wt_mouse_ratno$tspec_1 ~ 
                      expression_wt_mouse_ratno$tspec_2)$coefficients
  cor_value = cor.test(expression_wt_mouse_ratno$tspec_1, 
                       expression_wt_mouse_ratno$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  dev.off()
}
