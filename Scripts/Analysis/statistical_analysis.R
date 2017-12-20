'''
171220
Joaquim Claivaz

Statistical analysis
'''


####FUN####
classic_anco = function(species1, species2, considered_dataset){
  data_anco = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species1, '_', species2, considered_dataset))
  data_anco$status = as.character(data_anco$status)
  data_anco$status[data_anco$status != 'control'] = 'modif'
  data_anco$status = as.factor(data_anco$status)
  mod_anco = lm(data_anco$tspec_1 ~ data_anco$tspec_2 * data_anco$status)
  anova(mod_anco)
}

####ANCOVA analysis on the different datasets####
classic_anco(species1 = 'BOVIN', species2 = 'GORGO', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'BOVIN', species2 = 'HUMAN', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'BOVIN', species2 = 'MACMU', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'BOVIN', species2 = 'MONDO', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'BOVIN', species2 = 'MOUSE', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'BOVIN', species2 = 'PANTR', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'BOVIN', species2 = 'PIGXX', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'BOVIN', species2 = 'RATNO', considered_dataset = '_ortho_notfemale_nottestis_dataset')

classic_anco(species1 = 'GORGO', species2 = 'HUMAN', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'GORGO', species2 = 'MACMU', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'GORGO', species2 = 'MONDO', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'GORGO', species2 = 'MOUSE', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'GORGO', species2 = 'PANTR', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'GORGO', species2 = 'PIGXX', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'GORGO', species2 = 'RATNO', considered_dataset = '_ortho_notfemale_nottestis_dataset')

classic_anco(species2 = 'HUMAN', species1 = 'MACMU', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species2 = 'HUMAN', species1 = 'MONDO', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species2 = 'HUMAN', species1 = 'MOUSE', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species2 = 'HUMAN', species1 = 'PANTR', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species2 = 'HUMAN', species1 = 'PIGXX', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species2 = 'HUMAN', species1 = 'RATNO', considered_dataset = '_ortho_notfemale_nottestis_dataset')

classic_anco(species1 = 'MACMU', species2 = 'MONDO', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'MACMU', species2 = 'MOUSE', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'MACMU', species2 = 'PANTR', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'MACMU', species2 = 'PIGXX', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'MACMU', species2 = 'RATNO', considered_dataset = '_ortho_notfemale_nottestis_dataset')

classic_anco(species1 = 'MONDO', species2 = 'MOUSE', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'MONDO', species2 = 'PANTR', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'MONDO', species2 = 'PIGXX', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'MONDO', species2 = 'RATNO', considered_dataset = '_ortho_notfemale_nottestis_dataset')

classic_anco(species1 = 'MOUSE', species2 = 'PANTR', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'MOUSE', species2 = 'PIGXX', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'MOUSE', species2 = 'RATNO', considered_dataset = '_ortho_notfemale_nottestis_dataset')

classic_anco(species1 = 'PANTR', species2 = 'PIGXX', considered_dataset = '_ortho_notfemale_nottestis_dataset')
classic_anco(species1 = 'PANTR', species2 = 'RATNO', considered_dataset = '_ortho_notfemale_nottestis_dataset')

classic_anco(species1 = 'PIGXX', species2 = 'RATNO', considered_dataset = '_ortho_notfemale_nottestis_dataset')

####classic ancova - distribution and residual plot for human and mouse comparison####
species1 = 'HUMAN'
species2 = 'MOUSE'
considered_dataset = '_ortho_notfemale_dataset'
data_anco = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, considered_dataset))

data_anco$status = as.character(data_anco$status)
data_anco$status[data_anco$status != 'control'] = 'modif'
data_anco$status = as.factor(data_anco$status)
mod_anco = lm(data_anco$tspec_1 ~ data_anco$tspec_2 * data_anco$status)
qqnorm(residuals(mod_anco)); qqline(residuals(mod_anco))
plot(fitted.values(mod_anco), residuals(mod_anco), cex.main = 0.8, main = 'Residuals distribution in function of predicted values of the model (mod_anco)')
anova(mod_anco)

####distribution between control and modification group####
#MOUSE-HUMAN
species1 = 'HUMAN'
species2 = 'MOUSE'
considered_dataset = '_ortho_notfemale_dataset'
data_anco = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, considered_dataset))

hist(data_anco$tspec_1[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_1[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(data_anco$tspec_2[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_2[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

wilcox.test(data_anco$tspec_1[data_anco$status != 'control'],data_anco$tspec_1[data_anco$status == 'control'])
wilcox.test(data_anco$tspec_2[data_anco$status != 'control'],data_anco$tspec_2[data_anco$status == 'control'])

#MOUSE-RATNO
species1 = 'RATNO'
species2 = 'MOUSE'
considered_dataset = '_ortho_notfemale_dataset'
data_anco = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, considered_dataset))

hist(data_anco$tspec_1[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_1[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(data_anco$tspec_2[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_2[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

wilcox.test(data_anco$tspec_1[data_anco$status != 'control'],data_anco$tspec_1[data_anco$status == 'control'])
wilcox.test(data_anco$tspec_2[data_anco$status != 'control'],data_anco$tspec_2[data_anco$status == 'control'])

#HUMAN-MACMU
species1 = 'HUMAN'
species2 = 'MACMU'
considered_dataset = '_ortho_notfemale_dataset'
data_anco = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, considered_dataset))

hist(data_anco$tspec_1[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_1[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(data_anco$tspec_2[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_2[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

wilcox.test(data_anco$tspec_1[data_anco$status != 'control'],data_anco$tspec_1[data_anco$status == 'control'])
wilcox.test(data_anco$tspec_2[data_anco$status != 'control'],data_anco$tspec_2[data_anco$status == 'control'])

#HUMAN-BOVIN
species1 = 'HUMAN'
species2 = 'BOVIN'
considered_dataset = '_ortho_notfemale_dataset'
data_anco = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, considered_dataset))

hist(data_anco$tspec_1[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_1[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(data_anco$tspec_2[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_2[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

wilcox.test(data_anco$tspec_1[data_anco$status != 'control'],data_anco$tspec_1[data_anco$status == 'control'])
wilcox.test(data_anco$tspec_2[data_anco$status != 'control'],data_anco$tspec_2[data_anco$status == 'control'])

#RATNO-BOVIN
species1 = 'RATNO'
species2 = 'BOVIN'
considered_dataset = '_ortho_notfemale_dataset'
data_anco = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, considered_dataset))

hist(data_anco$tspec_1[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_1[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(data_anco$tspec_2[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), 
                   '\n Domain modification inferred from ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_2[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

wilcox.test(data_anco$tspec_1[data_anco$status != 'control'],data_anco$tspec_1[data_anco$status == 'control'])
wilcox.test(data_anco$tspec_2[data_anco$status != 'control'],data_anco$tspec_2[data_anco$status == 'control'])

####ANCOVA control and non parametric equivalent####
#Mann Whitney rank sum test (test if there's a difference in rank between modification group)
data_anco$status = as.character(data_anco$status)
data_anco$status[data_anco$status != 'control'] = 'modif'
hist(data_anco$tspec_1[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]), 
                   '\n in ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' in pairwise comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_1[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')
hist(data_anco$tspec_2[data_anco$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), 
                   '\n in ', gsub('[[:digit:]]','' ,data_anco$GeneID_1[1]),
                   ' and ', gsub('[[:digit:]]','' ,data_anco$GeneID_2[1]), ' in pairwise comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(data_anco$tspec_2[data_anco$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')
wilcox.test(data_anco$tspec_2 ~ data_anco$status)
wilcox.test(data_anco$tspec_1 ~ data_anco$status)

#paired t test (parametric and non parametric)
list_file = list.files('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/')
correlation_control = c()
correlation_modif = c()
j = 0
for (i in 1:length(list_file)){
  if (grepl('_ortho_notfemale_dataset', list_file[i])){
    j = j + 1
    data_tmp = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', list_file[i]))
    correlation_modif[j] = cor.test(data_tmp$tspec_2[data_tmp$status != 'control'], data_tmp$tspec_1[data_tmp$status != 'control'])$estimate[[1]]
    correlation_control[j] = cor.test(data_tmp$tspec_2[data_tmp$status == 'control'], data_tmp$tspec_1[data_tmp$status == 'control'])$estimate[[1]]
  }
}

correlation_data = data.frame(control = correlation_control, modif = correlation_modif)

hist(correlation_modif, breaks = 8, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of Pearson correlation factor'), xlab = 'R values', cex.main = 0.9, xlim = c(0,1), ylim = c(0,5))
hist(correlation_control, breaks = 8, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('topleft', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

t.test(correlation_data$control, correlation_data$modif)
wilcox.test(correlation_data$control, correlation_data$modif, paired = T)

####Gene length  bias control####
#previously the gene length (amino acid length) were inferred using EMBOSS (UNIX command)
#infoseq -auto -nocolumns -delimiter '\t' -only -noheading -name -length protein_sequence_BOVIN > protein_length
#infoseq -auto -nocolumns -delimiter '\t' -only -noheading -name -length protein_sequence_GORGO >> protein_length

##HUMAN & MOUSE
#load files
species1 = 'HUMAN'
species2 = 'MOUSE'
protein_length = read.table('D:/UNIL/Master/Master_Project/Data/OMA/protein_length')
data_genelength = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, '_ortho_notfemale_dataset'))
data_genelength = data_genelength[data_genelength$status == 'control',]

#data organization
goi_length = c(as.character(data_genelength$GeneID_1), as.character(data_genelength$GeneID_2))
goi_emboss_length = c()

for (i in 1:length(goi_length)){
  goi_emboss_length[i] = protein_length$V2[protein_length$V1 == goi_length[i]]
}

data_genelength$length_gene1 = goi_emboss_length[1:dim(data_genelength)[1]]
data_genelength$length_gene2 = goi_emboss_length[(1 + dim(data_genelength)[1]):(2 * dim(data_genelength)[1])]
data_genelength$dif_length = abs(data_genelength$length_gene1 - data_genelength$length_gene2)

data_genelength = data_genelength[order(data_genelength$dif_length),]
data_smalldiff = data_genelength[1:1000,]
data_bigdiff = data_genelength[(dim(data_genelength)[1] - 999):dim(data_genelength)[1],]

#graphical result
smoothScatter(data_smalldiff$tspec_1, data_smalldiff$tspec_2,
              xlab = '',
              ylab = '',
              main = 'Control group (1000 smallest gene length difference)', cex.main = 0.8, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,data_smalldiff$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,data_smalldiff$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(data_smalldiff$tspec_2 ~ data_smalldiff$tspec_1), col = 'red')
linear_param = lm(data_smalldiff$tspec_2 ~ data_smalldiff$tspec_1)$coefficients
cor_value = cor.test(data_smalldiff$tspec_2, data_smalldiff$tspec_1)$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(data_bigdiff$tspec_1, data_bigdiff$tspec_2,
              xlab = '',
              ylab = '',
              main = 'Control group (1000 biggest gene length difference)', cex.main = 0.8, cex.lab = 0.6)
title(xlab = paste0(gsub('[[:digit:]]','' ,data_bigdiff$GeneID_1[1]), ' Tspec values'), cex.lab =0.6, line = 2)
title(ylab = paste0(gsub('[[:digit:]]','' ,data_bigdiff$GeneID_2[1]), ' Tspec values'), cex.lab =0.6, line = 2)
abline(lm(data_bigdiff$tspec_2 ~ data_bigdiff$tspec_1), col = 'red')
linear_param = lm(data_bigdiff$tspec_2 ~ data_bigdiff$tspec_1)$coefficients
cor_value = cor.test(data_bigdiff$tspec_2, data_bigdiff$tspec_1)$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n r = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

####pie chart####
species1 = 'HUMAN'
species2 = 'MOUSE'
data_pie = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, '_ortho_notfemale_dataset'))

#keep only modified genes
data_pie = data_pie[data_pie$status != 'control',]
#keep only gene with shift from ubiquitous to specific
data_pie = data_pie[data_pie$shift %in% c('specific_to_ubiquitous', 'ubiquitous_to_specific'),]
#keep the tissue where the shift comes from/goes to
data_pie$tissueOI = NA
for (i in 1:dim(data_pie)[1]){
  if(data_pie$TspecF_1[i] == 'ubiquitous'){
    data_pie$tissueOI[i] = as.character(data_pie$spec_tissue_2[i])
  }else{
    data_pie$tissueOI[i] = as.character(data_pie$spec_tissue_1[i])
  }
}

#Visualization
labels_tmp = unique(data_pie$tissueOI)
percent_tmp = c()
for (i in 1:length(labels_tmp)){
  percent_tmp[i] = sum(data_pie$tissueOI == labels_tmp[i]) / length(data_pie$tissueOI)
}

pie(percent_tmp, labels = paste0(labels_tmp, ' ', round(percent_tmp, 3) * 100, '%'), main = paste0("Percentage of tissue where domain modification induces shifts between ubiquitous and specific\nDomain modification orthologs between ", gsub('[[:digit:]]','' ,data_pie$GeneID_1[1]), ' and ', gsub('[[:digit:]]','' ,data_pie$GeneID_2[1])), cex.main = 0.8, cex = 0.6)

#HUMAN-MACMU
species1 = 'HUMAN'
species2 = 'MACMU'
data_pie = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, '_ortho_notfemale_dataset'))

#keep only modified genes
data_pie = data_pie[data_pie$status != 'control',]
#keep only gene with shift from ubiquitous to specific
data_pie = data_pie[data_pie$shift %in% c('specific_to_ubiquitous', 'ubiquitous_to_specific'),]
#keep the tissue where the shift comes from/goes to
data_pie$tissueOI = NA
for (i in 1:dim(data_pie)[1]){
  if(data_pie$TspecF_1[i] == 'ubiquitous'){
    data_pie$tissueOI[i] = as.character(data_pie$spec_tissue_2[i])
  }else{
    data_pie$tissueOI[i] = as.character(data_pie$spec_tissue_1[i])
  }
}

#Visualization
labels_tmp = unique(data_pie$tissueOI)
percent_tmp = c()
for (i in 1:length(labels_tmp)){
  percent_tmp[i] = sum(data_pie$tissueOI == labels_tmp[i]) / length(data_pie$tissueOI)
}

pie(percent_tmp, labels = paste0(labels_tmp, ' ', round(percent_tmp, 3) * 100, '%'), main = paste0("Percentage of tissue where domain modification induces shifts between ubiquitous and specific\nDomain modification orthologs between ", gsub('[[:digit:]]','' ,data_pie$GeneID_1[1]), ' and ', gsub('[[:digit:]]','' ,data_pie$GeneID_2[1])), cex.main = 0.8, cex = 0.6)

#MOUSE-RATNO
species1 = 'RATNO'
species2 = 'MOUSE'
data_pie = read.csv(paste0('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/', species2, '_', species1, '_ortho_notfemale_dataset'))

#keep only modified genes
data_pie = data_pie[data_pie$status != 'control',]
#keep only gene with shift from ubiquitous to specific
data_pie = data_pie[data_pie$shift %in% c('specific_to_ubiquitous', 'ubiquitous_to_specific'),]
#keep the tissue where the shift comes from/goes to
data_pie$tissueOI = NA
for (i in 1:dim(data_pie)[1]){
  if(data_pie$TspecF_1[i] == 'ubiquitous'){
    data_pie$tissueOI[i] = as.character(data_pie$spec_tissue_2[i])
  }else{
    data_pie$tissueOI[i] = as.character(data_pie$spec_tissue_1[i])
  }
}

#Visualization
labels_tmp = unique(data_pie$tissueOI)
percent_tmp = c()
for (i in 1:length(labels_tmp)){
  percent_tmp[i] = sum(data_pie$tissueOI == labels_tmp[i]) / length(data_pie$tissueOI)
}

pie(percent_tmp, labels = paste0(labels_tmp, ' ', round(percent_tmp, 3) * 100, '%'), main = paste0("Percentage of tissue where domain modification induces shifts between ubiquitous and specific\nDomain modification orthologs between ", gsub('[[:digit:]]','' ,data_pie$GeneID_1[1]), ' and ', gsub('[[:digit:]]','' ,data_pie$GeneID_2[1])), cex.main = 0.8, cex = 0.6)

