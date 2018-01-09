'''
Joaquim Claivaz
171016

Tspec analysis plot results
'''
path_folder = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/'
different_species = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN')
central_species = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN')

regexp_list = c('_ortho_notfemale_dataset', '_ortho_onlymale_dataset',
                '_ortho_notfemale_onlybrain_dataset', '_ortho_onlymale_onlybrain_dataset',
                '_ortho_notfemale_nottestis_dataset', '_ortho_onlymale_nottestis_dataset',
                '_ortho_notfemale_onlybrain_nottestis_dataset', '_ortho_onlymale_onlybrain_nottestis_dataset',
                '_ortho_notfemale_onlypref_nottestis_dataset', '_ortho_onlymale_onlypref_nottestis_dataset',
                '_ortho_notfemale_onlypref_dataset', '_ortho_onlymale_onlypref_dataset') 

regexp_list = c('_para_notfemale_dataset', '_para_onlymale_dataset',
                '_para_notfemale_onlybrain_dataset', '_para_onlymale_onlybrain_dataset',
                '_para_notfemale_nottestis_dataset', '_para_onlymale_nottestis_dataset',
                '_para_notfemale_onlybrain_nottestis_dataset', '_para_onlymale_onlybrain_nottestis_dataset',
                '_para_notfemale_onlypref_nottestis_dataset', '_para_onlymale_onlypref_nottestis_dataset',
                '_para_notfemale_onlypref_dataset', '_para_onlymale_onlypref_dataset') 

####SEX factor consideration####
{
pdf(paste0('C:/Users/Claivaz/Desktop/', regexp_list[regexp_out] ,'.pdf'))
par(mfrow = c(5,4), mai=c(0.4,0.35,0.3,0.01))
specie_done = c()

species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                               '_ortho_notfemale_dataset'))
        
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

species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                               '_ortho_onlymale_dataset'))

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

#

species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'HUMAN',
                               '_ortho_notfemale_dataset'))

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

species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'HUMAN',
                               '_ortho_onlymale_dataset'))

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
#

species_data = read.csv(paste0(path_folder, 'GORGO', '_', 'MACMU',
                               '_ortho_notfemale_dataset'))

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

species_data = read.csv(paste0(path_folder, 'GORGO', '_', 'MACMU',
                               '_ortho_onlymale_dataset'))

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
#

species_data = read.csv(paste0(path_folder, 'MOUSE',
                               '_para_notfemale_dataset'))
species_data$X = NULL

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

species_data = read.csv(paste0(path_folder, 'MOUSE',
                               '_para_onlymale_dataset'))
species_data$X = NULL

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
#

species_data = read.csv(paste0(path_folder, 'BOVIN',
                               '_para_notfemale_dataset'))
species_data$X = NULL

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

species_data = read.csv(paste0(path_folder, 'BOVIN',
                               '_para_onlymale_dataset'))
species_data$X = NULL


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

dev.off()
}

####CEREBRAL tissue consideration####
{
pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
par(mfrow = c(4,6), mai=c(0.4,0.35,0.3,0.01))
specie_done = c()

species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                               '_ortho_notfemale_onlybrain_dataset'))

smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')


species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                               '_ortho_notfemale_onlypref_dataset'))

smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                               '_ortho_notfemale_dataset'))

smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')


species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'MACMU',
                               '_ortho_notfemale_onlybrain_dataset'))

smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')


species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'MACMU',
                               '_ortho_notfemale_onlypref_dataset'))

smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'MACMU',
                               '_ortho_notfemale_dataset'))

smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')


species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'HUMAN',
                               '_ortho_notfemale_onlybrain_dataset'))

smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')


species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'HUMAN',
                               '_ortho_notfemale_onlypref_dataset'))

smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'HUMAN',
                               '_ortho_notfemale_dataset'))

smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')


species_data = read.csv(paste0(path_folder, 'HUMAN',
                               '_para_notfemale_onlybrain_dataset'))
species_data$X = NULL


smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

species_data = read.csv(paste0(path_folder, 'HUMAN',
                               '_para_notfemale_onlypref_dataset'))
species_data$X = NULL


smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')



species_data = read.csv(paste0(path_folder, 'HUMAN',
                               '_para_notfemale_dataset'))
species_data$X = NULL


smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain modification group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = '',
              ylab = '',
              main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')


dev.off()
}
####TESTIS consideration####
{
pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
par(mfrow = c(5,4), mai=c(0.4,0.35,0.3,0.01))
specie_done = c()

species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                               '_ortho_notfemale_nottestis_dataset'))

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

species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                               '_ortho_notfemale_dataset'))

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
#

species_data = read.csv(paste0(path_folder, 'MONDO', '_', 'PANTR',
                               '_ortho_notfemale_nottestis_dataset'))

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

species_data = read.csv(paste0(path_folder, 'MONDO', '_', 'PANTR',
                               '_ortho_notfemale_dataset'))

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
#

species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'RATNO',
                               '_ortho_notfemale_nottestis_dataset'))

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

species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'RATNO',
                               '_ortho_notfemale_dataset'))

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
#

species_data = read.csv(paste0(path_folder, 'GORGO',
                               '_para_notfemale_nottestis_dataset'))
species_data$X = NULL


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

species_data = read.csv(paste0(path_folder, 'GORGO',
                               '_para_notfemale_dataset'))
species_data$X = NULL


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
#

species_data = read.csv(paste0(path_folder, 'HUMAN',
                               '_para_notfemale_nottestis_dataset'))
species_data$X = NULL


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

species_data = read.csv(paste0(path_folder, 'HUMAN',
                               '_para_notfemale_dataset'))
species_data$X = NULL


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

dev.off()
}

####Paralog sorting####
{
pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
par(mfrow = c(4,3), mai=c(0.4,0.35,0.3,0.01))


species_data = read.csv(paste0(path_folder, 'HUMAN',
                               '_para_notfemale_dataset'))
species_data$X = NULL


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


species_data = read.csv(paste0(path_folder, 'HUMAN',
                               '_para_notfemale_dataset'))
species_data$X = NULL


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



species_data = read.csv(paste0(path_folder, 'MOUSE',
                               '_para_notfemale_dataset'))
species_data$X = NULL


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


species_data = read.csv(paste0(path_folder, 'MOUSE',
                               '_para_notfemale_dataset'))
species_data$X = NULL

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


species_data = read.csv(paste0(path_folder, 'BOVIN',
                               '_para_notfemale_dataset'))
species_data$X = NULL


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


species_data = read.csv(paste0(path_folder, 'BOVIN',
                               '_para_notfemale_dataset'))
species_data$X = NULL


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


species_data = read.csv(paste0(path_folder, 'MACMU',
                               '_para_notfemale_dataset'))
species_data$X = NULL


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


species_data = read.csv(paste0(path_folder, 'MACMU',
                               '_para_notfemale_dataset'))
species_data$X = NULL

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


dev.off()
}

####TISSUE_SELECTION####
tissue_selected = 'notfemale_dataset'
tissue_selected = 'notfemale_nottestis_dataset'

####HOMINIDAE#####
##ortholog
{
pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,4), mai=c(0.4,0.35,0.3,0.01))
  
  
species_data = read.csv(paste0(path_folder, 'GORGO', '_', 'HUMAN',
                               '_ortho_', tissue_selected))

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

species_data = read.csv(paste0(path_folder, 'GORGO', '_', 'MACMU',
                               '_ortho_', tissue_selected))

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

species_data = read.csv(paste0(path_folder, 'GORGO', '_', 'PANTR',
                               '_ortho_', tissue_selected))

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

species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'HUMAN',
                               '_ortho_', tissue_selected))

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

species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'PANTR',
                               '_ortho_', tissue_selected))

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
species_data = read.csv(paste0(path_folder, 'PANTR', '_', 'HUMAN',
                               '_ortho_', tissue_selected))

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

dev.off()
}


{
pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
par(mfrow = c(6,4), mai=c(0.4,0.35,0.3,0.01))


species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'RATNO',
                               '_ortho_', tissue_selected))

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

dev.off()
}

####HOMINIDAE MURIDAE####
{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,4), mai=c(0.4,0.35,0.3,0.01))
  
  
  species_data = read.csv(paste0(path_folder, 'GORGO', '_', 'MOUSE',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'GORGO', '_', 'RATNO',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'MOUSE',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'RATNO',
                                 '_ortho_', tissue_selected))
  
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
  
  
  species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'PANTR',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'PANTR', '_', 'RATNO',
                                 '_ortho_', tissue_selected))
  
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
  
  
  
  species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'HUMAN',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'RATNO', '_', 'HUMAN',
                                 '_ortho_', tissue_selected))
  
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
  dev.off()
}

species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'HUMAN',
                               '_ortho_', tissue_selected))

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

{
pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,3), mai=c(0.4,0.35,0.3,0.01))
  
species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'RATNO',
                               '_ortho_', tissue_selected))

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
dev.off()
}


#####HUMAN-MOUSE-other####

{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,4), mai=c(0.4,0.35,0.3,0.01))
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'MOUSE',
                                 '_ortho_', tissue_selected))
  
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
  
  
  
  species_data = read.csv(paste0(path_folder, 'MONDO', '_', 'HUMAN',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'MONDO', '_', 'MOUSE',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'PIGXX', '_', 'HUMAN',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'MOuSE', '_', 'PIGXX',
                                 '_ortho_', tissue_selected))
  
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
  
  
  
  dev.off()
}

#####other####

{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,4), mai=c(0.4,0.35,0.3,0.01))
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'GORGO',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'MACMU',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'PANTR',
                                 '_ortho_', tissue_selected))
  
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
 
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'RATNO',
                                 '_ortho_', tissue_selected))
  
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
  
  
  species_data = read.csv(paste0(path_folder, 'GORGO', '_', 'MONDO',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'MONDO',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'MONDO', '_', 'PANTR',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'MONDO', '_', 'RATNO',
                                 '_ortho_', tissue_selected))
  
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
  
  
  species_data = read.csv(paste0(path_folder, 'GORGO', '_', 'PIGXX',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'PIGXX',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'PANTR', '_', 'PIGXX',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'PIGXX', '_', 'RATNO',
                                 '_ortho_', tissue_selected))
  
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
  
  
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'MONDO',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'PIGXX',
                                 '_ortho_', tissue_selected))
  
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
  
  species_data = read.csv(paste0(path_folder, 'MONDO', '_', 'PIGXX',
                                 '_ortho_', tissue_selected))
  
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
  
  
  dev.off()
}

#####PARALOG HOMINIDAE####
{
pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
par(mfrow = c(6,2), mai=c(0.4,0.35,0.3,0.01))
  

species_data = read.csv(paste0(path_folder, 'MACMU',
                               '_para_', tissue_selected))
species_data$X = NULL


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




species_data = read.csv(paste0(path_folder, 'GORGO',
                               '_para_', tissue_selected))
species_data$X = NULL

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


species_data = read.csv(paste0(path_folder, 'PANTR',
                               '_para_', tissue_selected))
species_data$X = NULL


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

species_data = read.csv(paste0(path_folder, 'HUMAN',
                               '_para_', tissue_selected))
species_data$X = NULL


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

dev.off()
}

####PARALOG muridae####
{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,2), mai=c(0.4,0.35,0.3,0.01))
  
  
  species_data = read.csv(paste0(path_folder, 'MOUSE',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  
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
  
  
  
  
  species_data = read.csv(paste0(path_folder, 'RATNO',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  
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
  
  dev.off()
}  



#####PARALOG other####
{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,2), mai=c(0.4,0.35,0.3,0.01))
  
  
  species_data = read.csv(paste0(path_folder, 'BOVIN',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  
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
  
  
  
  
  species_data = read.csv(paste0(path_folder, 'MONDO',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  
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
  
  
  species_data = read.csv(paste0(path_folder, 'PIGXX',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  
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
  
  
  
  
  
  dev.off()
}  

####ORTHOLOG pairs effect loss####
{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,4), mai=c(0.4,0.35,0.3,0.01))
  
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                                 '_ortho_', tissue_selected))
  

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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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

species_data = read.csv(paste0(path_folder, 'MACMU', '_', 'HUMAN',
                               '_ortho_', tissue_selected))


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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'HUMAN',
                               '_ortho_', tissue_selected))


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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
species_data = read.csv(paste0(path_folder, 'RATNO', '_', 'HUMAN',
                               '_ortho_', tissue_selected))


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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'MOUSE',
                               '_ortho_', tissue_selected))


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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'RATNO',
                               '_ortho_', tissue_selected))


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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

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

dev.off()
}  

####PARALOG pairs effect loss ref gene####
{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,4), mai=c(0.4,0.35,0.3,0.01))
  
  species_data = read.csv(paste0(path_folder, 'BOVIN',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'GORGO',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'HUMAN',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'MACMU',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'MONDO',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'MOUSE',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'PANTR',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'PIGXX',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'RATNO',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5, xlim = c(0,1), ylim = c(0,1))
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, xlim = c(0,1), ylim = c(0,1), cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  
  dev.off()
}  


####PARALOG pairs effect loss all modif####
{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,4), mai=c(0.4,0.35,0.3,0.01))
  
  species_data = read.csv(paste0(path_folder, 'BOVIN',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'GORGO',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'HUMAN',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'MACMU',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'MONDO',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'MOUSE',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'PANTR',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'PIGXX',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  species_data = read.csv(paste0(path_folder, 'RATNO',
                                 '_para_', tissue_selected))
  species_data$X = NULL
  
  smoothScatter(species_data$tspec_1[species_data$status == 'f-1'],
                species_data$tspec_2[species_data$status == 'f-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: f-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'f-1'] ~
              species_data$tspec_1[species_data$status == 'f-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'f-1'] ~
                      species_data$tspec_1[species_data$status == 'f-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'f-1'],
                       species_data$tspec_1[species_data$status == 'f-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'int-1'],
                species_data$tspec_2[species_data$status == 'int-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: int-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'int-1'] ~
              species_data$tspec_1[species_data$status == 'int-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'int-1'] ~
                      species_data$tspec_1[species_data$status == 'int-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'int-1'],
                       species_data$tspec_1[species_data$status == 'int-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'b-1'],
                species_data$tspec_2[species_data$status == 'b-1'],
                xlab = '',
                ylab = '',
                main = paste0('Domain modification group: b-1'), cex.main = 0.5)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'b-1'] ~
              species_data$tspec_1[species_data$status == 'b-1']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'b-1'] ~
                      species_data$tspec_1[species_data$status == 'b-1'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'b-1'],
                       species_data$tspec_1[species_data$status == 'b-1'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                xlab = '',
                ylab = '',
                main = 'Domain control group', cex.main = 0.5, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' ref Tspec values'), cex.lab =0.6, line = 2)
  title(ylab = paste0('Other ' ,gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
  abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
  linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
  cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  
  dev.off()
}  


####Paralog- inter species ####
path_folder = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/'
different_species = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN')
central_species = c('BOVIN', 'GORGO', 'MACMU', 'MONDO', 'MOUSE', 'PANTR', 'PIGXX', 'RATNO', 'HUMAN')
regexp_list = c('_paraortho_notfemale_withouttestis_dataset', '_paraortho_notfemale_withtestis_dataset') 

for (regexp_out in 1:length(regexp_list)){
  pdf(paste0('C:/Users/Claivaz/Desktop/tmp_cor', regexp_list[regexp_out] ,'.pdf'))
  par(mfrow = c(6,4), mai=c(0.4,0.35,0.3,0.01))
  specie_done = c()
  
  for (sp1 in 1:length(central_species)){
    specie_done = c(specie_done, central_species[sp1])
    
    for (sp2 in 1:length(different_species)){
      
      if (central_species[sp1] != different_species[sp2] & !(different_species[sp2] %in% specie_done)){
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
        
        smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
                        xlab = '',
                        ylab = '',
                        main = 'Domain control group', cex.main = 0.8, cex.lab = 0.6)
          title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
          options(show.error.messages = FALSE)
          if(class(try(abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ 
                                 species_data$tspec_1[species_data$status == 'control']), col = 'red'))) != 'try-error'){
            abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ 
                        species_data$tspec_1[species_data$status == 'control']), col = 'red')
            if(class(try(cor.test(species_data$tspec_2[species_data$status == 'control'], 
                                  species_data$tspec_1[species_data$status == 'control'])$estimate[[1]])) != 'try-error'){
              linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
              cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
              text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
            }
          }
          options(show.error.messages = TRUE)
          
          
        
      }
    }
  }
  dev.off()
}


####longest domain group orthologs comparisons####
{
  pdf(paste0('C:/Users/Claivaz/Desktop/', 'tmp' ,'.pdf'))
  par(mfrow = c(6,3), mai=c(0.4,0.35,0.3,0.01))
  
  
  species_data = read.csv(paste0(path_folder, 'PANTR', '_', 'HUMAN',
                                 '_ortho_', tissue_selected))
  longer_sp1 = species_data[species_data$len_1 > species_data$len_2,]
  longer_sp2 = species_data[species_data$len_1 < species_data$len_2,]
  
  smoothScatter(longer_sp1$tspec_1, 
                longer_sp1$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp1$tspec_2 ~ 
              longer_sp1$tspec_1), col = 'red')
  linear_param = lm(longer_sp1$tspec_2 ~ 
                      longer_sp1$tspec_1)$coefficients
  cor_value = cor.test(longer_sp1$tspec_1,
                       longer_sp1$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(longer_sp2$tspec_1, 
                longer_sp2$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp2$tspec_2 ~ 
              longer_sp2$tspec_1), col = 'red')
  linear_param = lm(longer_sp2$tspec_2 ~ 
                      longer_sp2$tspec_1)$coefficients
  cor_value = cor.test(longer_sp2$tspec_2,
                       longer_sp2$tspec_1)$estimate[[1]]
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
  
  species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'HUMAN',
                                 '_ortho_', tissue_selected))
  longer_sp1 = species_data[species_data$len_1 > species_data$len_2,]
  longer_sp2 = species_data[species_data$len_1 < species_data$len_2,]
  
  smoothScatter(longer_sp1$tspec_1, 
                longer_sp1$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp1$tspec_2 ~ 
              longer_sp1$tspec_1), col = 'red')
  linear_param = lm(longer_sp1$tspec_2 ~ 
                      longer_sp1$tspec_1)$coefficients
  cor_value = cor.test(longer_sp1$tspec_1,
                       longer_sp1$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(longer_sp2$tspec_1, 
                longer_sp2$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp2$tspec_2 ~ 
              longer_sp2$tspec_1), col = 'red')
  linear_param = lm(longer_sp2$tspec_2 ~ 
                      longer_sp2$tspec_1)$coefficients
  cor_value = cor.test(longer_sp2$tspec_2,
                       longer_sp2$tspec_1)$estimate[[1]]
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
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'HUMAN',
                                 '_ortho_', tissue_selected))
  longer_sp1 = species_data[species_data$len_1 > species_data$len_2,]
  longer_sp2 = species_data[species_data$len_1 < species_data$len_2,]
  
  smoothScatter(longer_sp1$tspec_1, 
                longer_sp1$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp1$tspec_2 ~ 
              longer_sp1$tspec_1), col = 'red')
  linear_param = lm(longer_sp1$tspec_2 ~ 
                      longer_sp1$tspec_1)$coefficients
  cor_value = cor.test(longer_sp1$tspec_1,
                       longer_sp1$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(longer_sp2$tspec_1, 
                longer_sp2$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp2$tspec_2 ~ 
              longer_sp2$tspec_1), col = 'red')
  linear_param = lm(longer_sp2$tspec_2 ~ 
                      longer_sp2$tspec_1)$coefficients
  cor_value = cor.test(longer_sp2$tspec_2,
                       longer_sp2$tspec_1)$estimate[[1]]
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
  
  species_data = read.csv(paste0(path_folder, 'MOUSE', '_', 'RATNO',
                                 '_ortho_', tissue_selected))
  longer_sp1 = species_data[species_data$len_1 > species_data$len_2,]
  longer_sp2 = species_data[species_data$len_1 < species_data$len_2,]
  
  smoothScatter(longer_sp1$tspec_1, 
                longer_sp1$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp1$tspec_2 ~ 
              longer_sp1$tspec_1), col = 'red')
  linear_param = lm(longer_sp1$tspec_2 ~ 
                      longer_sp1$tspec_1)$coefficients
  cor_value = cor.test(longer_sp1$tspec_1,
                       longer_sp1$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(longer_sp2$tspec_1, 
                longer_sp2$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp2$tspec_2 ~ 
              longer_sp2$tspec_1), col = 'red')
  linear_param = lm(longer_sp2$tspec_2 ~ 
                      longer_sp2$tspec_1)$coefficients
  cor_value = cor.test(longer_sp2$tspec_2,
                       longer_sp2$tspec_1)$estimate[[1]]
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
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'MOUSE',
                                 '_ortho_', tissue_selected))
  longer_sp1 = species_data[species_data$len_1 > species_data$len_2,]
  longer_sp2 = species_data[species_data$len_1 < species_data$len_2,]
  
  smoothScatter(longer_sp1$tspec_1, 
                longer_sp1$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp1$tspec_2 ~ 
              longer_sp1$tspec_1), col = 'red')
  linear_param = lm(longer_sp1$tspec_2 ~ 
                      longer_sp1$tspec_1)$coefficients
  cor_value = cor.test(longer_sp1$tspec_1,
                       longer_sp1$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(longer_sp2$tspec_1, 
                longer_sp2$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp2$tspec_2 ~ 
              longer_sp2$tspec_1), col = 'red')
  linear_param = lm(longer_sp2$tspec_2 ~ 
                      longer_sp2$tspec_1)$coefficients
  cor_value = cor.test(longer_sp2$tspec_2,
                       longer_sp2$tspec_1)$estimate[[1]]
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
  
  species_data = read.csv(paste0(path_folder, 'BOVIN', '_', 'RATNO',
                                 '_ortho_', tissue_selected))
  longer_sp1 = species_data[species_data$len_1 > species_data$len_2,]
  longer_sp2 = species_data[species_data$len_1 < species_data$len_2,]
  
  smoothScatter(longer_sp1$tspec_1, 
                longer_sp1$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp1$tspec_2 ~ 
              longer_sp1$tspec_1), col = 'red')
  linear_param = lm(longer_sp1$tspec_2 ~ 
                      longer_sp1$tspec_1)$coefficients
  cor_value = cor.test(longer_sp1$tspec_1,
                       longer_sp1$tspec_2)$estimate[[1]]
  text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')
  
  smoothScatter(longer_sp2$tspec_1, 
                longer_sp2$tspec_2,
                xlab = '',
                ylab = '',
                main = 'Domain modification group', cex.main = 0.8, cex.lab = 0.6)
  title(xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values (shortest pair)'), cex.lab =0.6, line = 2)
  title(ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values (longest pair)'), cex.lab =0.6, line = 2)
  abline(lm(longer_sp2$tspec_2 ~ 
              longer_sp2$tspec_1), col = 'red')
  linear_param = lm(longer_sp2$tspec_2 ~ 
                      longer_sp2$tspec_1)$coefficients
  cor_value = cor.test(longer_sp2$tspec_2,
                       longer_sp2$tspec_1)$estimate[[1]]
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
  
  dev.off()
}
