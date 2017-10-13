'''
Joaquim Claivaz
171003, last modification 171013

Tspec analysis in mammals
'''
####FUN####
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
####

#####Tspec analysis#####
species_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/BOVIN_GORGO_ortho_notfemale_dataset')
species_data$X = NULL
####Ortholog####
##status
#Effect of domain modification (status) on Tspec values 
hist(species_data$tspec_1[species_data$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), '\n', gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(species_data$tspec_1[species_data$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')

hist(species_data$tspec_2[species_data$status != 'control'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), '\n', gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(species_data$tspec_2[species_data$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')

#H0: no differences amongst distribution of tspec values in function of domain modification status
ks.test(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_1[species_data$status == 'control'])
ks.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_2[species_data$status == 'control'])
#

#Effect of domain modification (status) on specificity factor (ubiquitous / specificity)
'''
c_s = sum(species_data$TspecF_1[species_data$status == 'control'] == 'specific')
c_u = sum(species_data$TspecF_1[species_data$status == 'control'] == 'ubiquitous')
m_s = sum(species_data$TspecF_1[species_data$status != 'control'] == 'specific')
m_u = sum(species_data$TspecF_1[species_data$status != 'control'] == 'ubiquitous')

barplot(c(c(c_s / (c_s + c_u), c_u / (c_s + c_u)),
          c(m_s / (m_s + m_u), m_u / (m_s + m_u))),
        names.arg = rep(c('specific', 'ubiquitous'), 2),
        cex.names = 0.6, col = rep(c('blue', 'red'), each = 2), beside = T,
        main = paste0('Pair proportion in each group of specificity factor\n', gsub('[[:digit:]]', '', species_data$GeneID_1[1])), ylab = 'pair proportion')
legend('bottomright', c( "domain modification", "control"), fill = c('red','blue'),
       cex = 0.6, horiz = F)

c_s = sum(species_data$TspecF_2[species_data$status == 'control'] == 'specific')
c_u = sum(species_data$TspecF_2[species_data$status == 'control'] == 'ubiquitous')
m_s = sum(species_data$TspecF_2[species_data$status != 'control'] == 'specific')
m_u = sum(species_data$TspecF_2[species_data$status != 'control'] == 'ubiquitous')

barplot(c(c(c_s / (c_s + c_u), c_u / (c_s + c_u)),
          c(m_s / (m_s + m_u), m_u / (m_s + m_u))),
        names.arg = rep(c('specific', 'ubiquitous'), 2),
        cex.names = 0.6, col = rep(c('blue', 'red'), each = 2), beside = T,
        main = paste0('Pair proportion in each group of specificity factor\n', gsub('[[:digit:]]', '', species_data$GeneID_2[1])), ylab = 'pair proportion')
legend('bottomright', c( "domain modification", "control"), fill = c('red','blue'),
       cex = 0.6, horiz = F)


#HO: no difference amongst specific and ubiquitous factor proportion in function of domain modification status
chisq.test(species_data$TspecF_1, species_data$status, correct = F)
chisq.test(species_data$TspecF_2, species_data$status, correct = F)
'''
#

#Correlation between tspec in function of status
smoothScatter(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
               xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'),
               ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'),
               main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,species_data$GeneID_2[1]),
                             ' and ', gsub('[[:digit:]]','' ,species_data$GeneID_1[1]),
                             '\nDomain modification group'), cex.main = 0.8)
abline(lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status != 'control'] ~ species_data$tspec_1[species_data$status != 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status != 'control'], species_data$tspec_1[species_data$status != 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n R2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
               xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'),
               ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'),
               main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,species_data$GeneID_2[1]),
                             ' and ', gsub('[[:digit:]]','' ,species_data$GeneID_1[1]),
                             '\nDomain control group'), cex.main = 0.8)
abline(lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control']), col = 'red')
linear_param = lm(species_data$tspec_2[species_data$status == 'control'] ~ species_data$tspec_1[species_data$status == 'control'])$coefficients
cor_value = cor.test(species_data$tspec_2[species_data$status == 'control'], species_data$tspec_1[species_data$status == 'control'])$estimate[[1]]
text(0.7, 0.2, paste0(round(linear_param[1], digits = 3), ' + ', round(linear_param[2], digits = 3), ' x = y\n R2 = ', round(cor_value, digits = 3)), cex = 0.8, col = 'red')


#model_ancova = lm(species_data$tspec_1 ~ species_data$tspec_2 + species_data$status)
#qqnorm(residuals(model_ancova)); qqline(residuals(model_ancova))
#plot(fitted.values(model_ancova), residuals(model_ancova))
#summary(model_ancova)

cor.test(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'])
cor.test(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'])

#H0: no difference in the correlation factor in tspec values pairwise comparison in function of the considered group
cor.diff.test(species_data$tspec_1[species_data$status != 'control'], species_data$tspec_2[species_data$status != 'control'],
              species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'])
#

#study proportion of specific factor in function of status
'''
barplot(c(table(species_data$spec_tissue_1[species_data$status != 'control']) / sum(table(species_data$spec_tissue_1[species_data$status != 'control'])),
          table(species_data$spec_tissue_1[species_data$status == 'control']) / sum(table(species_data$spec_tissue_1[species_data$status == 'control']))),
        cex.names = 0.6, las = 2,
        col =  rep(c('red', 'blue'), each =length(levels(species_data$spec_tissue_1))),
        main = paste0('Pair proportion in each group of specificity factor\nIn ', gsub('[[:digit:]]', '', species_data$GeneID_1[1])), ylab = 'pair proportion')
legend('topleft', c( "domain modification", "control"), fill = c('red','blue'),
       cex = 0.6, horiz = F)

barplot(c(table(species_data$spec_tissue_2[species_data$status != 'control']) / sum(table(species_data$spec_tissue_2[species_data$status != 'control'])),
          table(species_data$spec_tissue_2[species_data$status == 'control']) / sum(table(species_data$spec_tissue_2[species_data$status == 'control']))),
          cex.names = 0.6, col = rep(c('red', 'blue'), each =length(levels(species_data$spec_tissue_2))), las = 2,
          main = paste0('Pair proportion in each group of specificity factor\nIn ', gsub('[[:digit:]]', '', species_data$GeneID_2[1])), ylab = 'pair proportion')
legend('topleft', c( "domain modification", "control"), fill = c('red','blue'),
       cex = 0.6, horiz = F)

#H0: no difference amongst the proportion of specificity factor in function of domain modification
chisq.test(species_data$spec_tissue_1, species_data$status, correct = F)
chisq.test(species_data$spec_tissue_2, species_data$status, correct = F)
'''
#

#study shift in specific factor in function of status
'''
barplot(c(table(species_data$shift[species_data$status != 'control']) / sum(table(species_data$shift[species_data$status != 'control'])),
          table(species_data$shift[species_data$status == 'control']) / sum(table(species_data$shift[species_data$status == 'control']))),
        cex.names = 0.6, col = rep(c('red', 'blue'), each = 5), las = 2, beside = T,
        main = paste0('Pair proportion in each group of shift specificity factor\n', gsub('[[:digit:]]', '', species_data$GeneID_1[1]), ' vs ', gsub('[[:digit:]]', '', species_data$GeneID_2[1])), ylab = 'pair proportion')
legend('top', c( "domain modification", "control"), fill = c('red','blue'),
       cex = 0.6, horiz = F)
'''

pie(table(species_data$shift[species_data$status != 'control']) / sum(table(species_data$shift[species_data$status != 'control'])),
        cex.names = 0.6 , 
        main = paste0('Pair proportion in modification group of shift specificity factor\n', gsub('[[:digit:]]', '', species_data$GeneID_1[1]), ' vs ', gsub('[[:digit:]]', '', species_data$GeneID_2[1])))

pie(table(species_data$shift[species_data$status == 'control']) / sum(table(species_data$shift[species_data$status == 'control'])),
    cex.names = 0.6 , 
    main = paste0('Pair proportion in control group of shift specificity factor\n', gsub('[[:digit:]]', '', species_data$GeneID_1[1]), ' vs ', gsub('[[:digit:]]', '', species_data$GeneID_2[1])))

#H0: no difference in proportion among control and modification
chisq.test(species_data$shift, species_data$status, correct = F)

#H0: no difference in proportion of ortholog evolving in specificity ('neo- / subfunctionalization) in function of domain modification
#Shift specificity: specificity_shift, ubiquitous_to_specific
#No shift specificity: specific_no_shift, specific_to_ubiquitous, ubiquitous_no_shift
chisq.test(c(rep(1, length(species_data$shift[species_data$shift == 'ubiquitous_to_specific' | species_data$shift == 'specificity_shift'])),
             rep(0, length(species_data$shift[!(species_data$shift == 'ubiquitous_to_specific' | species_data$shift == 'specificity_shift')]))),
           c(as.factor(species_data$status[species_data$shift == 'ubiquitous_to_specific' | species_data$shift == 'specificity_shift']),
             as.factor(species_data$status[!(species_data$shift == 'ubiquitous_to_specific' | species_data$shift == 'specificity_shift')])), correct = F)

#Shift specificity: specificity_shift, ubiquitous_to_specific, specific_to_ubiquitous
#No shift specificity: specific_no_shift, ubiquitous_no_shift
chisq.test(c(rep(1, length(species_data$shift[species_data$shift == 'ubiquitous_to_specific' | species_data$shift == 'specificity_shift' | species_data$shift == 'specific_to_ubiquitous'])), 
             rep(0, length(species_data$shift[!(species_data$shift == 'ubiquitous_to_specific' | species_data$shift == 'specificity_shift' | species_data$shift == 'specific_to_ubiquitous')]))), 
           c(as.factor(species_data$status[species_data$shift == 'ubiquitous_to_specific' | species_data$shift == 'specificity_shift' | species_data$shift == 'specific_to_ubiquitous']), 
             as.factor(species_data$status[!(species_data$shift == 'ubiquitous_to_specific' | species_data$shift == 'specificity_shift' | species_data$shift == 'specific_to_ubiquitous')])), correct = F)
##
####

###Effect of longer domain on the different estimators####
#on Tspec values 
hist(species_data$tspec_1[species_data$longer_domain_specie == 'specie1'],
     breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5),
     main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n',
                   gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' and ',gsub('[[:digit:]]','' ,species_data$GeneID_2[1]),
                   ' comparison\nlonger domain in ', gsub('[[:digit:]]', '', species_data$GeneID_1[1])),
     xlab = 'Tspec value', cex.main = 0.8)
hist(species_data$tspec_2[species_data$longer_domain_specie == 'specie1'], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c(paste0("tspec in ", gsub('[[:digit:]]', '', species_data$GeneID_1[1])), paste0("tpsec in ", gsub('[[:digit:]]', '', species_data$GeneID_2[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(species_data$tspec_1[species_data$longer_domain_specie == 'specie2'],
     breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5),
      main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n',
                   gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' and ',gsub('[[:digit:]]','' ,species_data$GeneID_2[1]),
                   ' comparison\nlonger domain in ', gsub('[[:digit:]]', '', species_data$GeneID_2[1])),
     xlab = 'Tspec value', cex.main = 0.8)
hist(species_data$tspec_2[species_data$longer_domain_specie == 'specie2'], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c(paste0("tspec in ", gsub('[[:digit:]]', '', species_data$GeneID_1[1])), paste0("tpsec in ", gsub('[[:digit:]]', '', species_data$GeneID_2[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(species_data$tspec_1[species_data$longer_domain_specie == 'specie1'],
     breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5),
     main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n',
                   gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' and ',gsub('[[:digit:]]','' ,species_data$GeneID_2[1]),
                   ' comparison\ntspec in ', gsub('[[:digit:]]', '', species_data$GeneID_1[1])),
     xlab = 'Tspec value', cex.main = 0.8)
hist(species_data$tspec_1[species_data$longer_domain_specie == 'specie2'], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c(paste0("longer domain in ", gsub('[[:digit:]]', '', species_data$GeneID_1[1])), paste0("longer domain in ", gsub('[[:digit:]]', '', species_data$GeneID_2[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')

hist(species_data$tspec_2[species_data$longer_domain_specie == 'specie1'],
     breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5),
     main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n',
                   gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' and ',gsub('[[:digit:]]','' ,species_data$GeneID_2[1]),
                   ' comparison\ntspec in ', gsub('[[:digit:]]', '', species_data$GeneID_2[1])),
     xlab = 'Tspec value', cex.main = 0.8)
hist(species_data$tspec_2[species_data$longer_domain_specie == 'specie2'], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c(paste0("longer domain in ", gsub('[[:digit:]]', '', species_data$GeneID_1[1])), paste0("longer domain in ", gsub('[[:digit:]]', '', species_data$GeneID_2[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')

#H0: no difference in the tspec values distribution in function of longer domain
ks.test(species_data$tspec_1[species_data$longer_domain_specie == 'specie2' & !(is.na(species_data$tspec_1[species_data$longer_domain_specie == 'specie2']))],
        species_data$tspec_1[species_data$longer_domain_specie == 'specie1' & !(is.na(species_data$tspec_1[species_data$longer_domain_specie == 'specie1']))])
ks.test(species_data$tspec_2[species_data$longer_domain_specie == 'specie2' & !(is.na(species_data$tspec_2[species_data$longer_domain_specie == 'specie2']))],
        species_data$tspec_2[species_data$longer_domain_specie == 'specie1' & !(is.na(species_data$tspec_2[species_data$longer_domain_specie == 'specie1']))])

#H0: no difference in tspec value distribution between longer domain and control group
ks.test(species_data$tspec_1[species_data$longer_domain_specie == 'specie2'],
        species_data$tspec_1[species_data$status == 'control'])
ks.test(species_data$tspec_1[species_data$longer_domain_specie == 'specie1'],
        species_data$tspec_1[species_data$status == 'control'])
ks.test(species_data$tspec_2[species_data$longer_domain_specie == 'specie2'],
        species_data$tspec_2[species_data$status == 'control'])
ks.test(species_data$tspec_2[species_data$longer_domain_specie == 'specie1'],
        species_data$tspec_2[species_data$status == 'control'])

smoothScatter(species_data$tspec_1[species_data$status != 'control' & species_data$longer_domain_specie == 'specie1'], 
              species_data$tspec_2[species_data$status != 'control'  & species_data$longer_domain_specie == 'specie1'],
              xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'),
              ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'),
              main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,species_data$GeneID_2[1]),
                            ' and ', gsub('[[:digit:]]','' ,species_data$GeneID_1[1]),
                            '\nDomain modification group, longer domain in ', gsub('[[:digit:]]', '', species_data$GeneID_1[1])), cex.main = 0.8)

smoothScatter(species_data$tspec_1[species_data$status != 'control' & species_data$longer_domain_specie == 'specie2'], 
              species_data$tspec_2[species_data$status != 'control'  & species_data$longer_domain_specie == 'specie2'],
              xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'),
              ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'),
              main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,species_data$GeneID_2[1]),
                            ' and ', gsub('[[:digit:]]','' ,species_data$GeneID_1[1]),
                            '\nDomain modification group, longer domain in ', gsub('[[:digit:]]', '', species_data$GeneID_2[1])), cex.main = 0.8)

smoothScatter(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'],
              xlab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_1[1]), ' Tspec values'),
              ylab = paste0(gsub('[[:digit:]]','' ,species_data$GeneID_2[1]), ' Tspec values'),
              main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,species_data$GeneID_2[1]),
                            ' and ', gsub('[[:digit:]]','' ,species_data$GeneID_1[1]),
                            '\nDomain control group'), cex.main = 0.8)

cor.test(species_data$tspec_1[species_data$status != 'control' & species_data$longer_domain_specie == 'specie2'], species_data$tspec_2[species_data$status != 'control'& species_data$longer_domain_specie == 'specie2'])
cor.test(species_data$tspec_1[species_data$status != 'control' & species_data$longer_domain_specie == 'specie1'], species_data$tspec_2[species_data$status != 'control'& species_data$longer_domain_specie == 'specie1'])
cor.test(species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'])

#H0: no difference in the correlation factor in tspec values pairwise comparison in function of the considered group
cor.diff.test(species_data$tspec_1[species_data$status != 'control' & species_data$longer_domain_specie == 'specie2'], species_data$tspec_2[species_data$status != 'control' & species_data$longer_domain_specie == 'specie2'],
              species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'])
cor.diff.test(species_data$tspec_1[species_data$status != 'control' & species_data$longer_domain_specie == 'specie1'], species_data$tspec_2[species_data$status != 'control' & species_data$longer_domain_specie == 'specie1'],
              species_data$tspec_1[species_data$status == 'control'], species_data$tspec_2[species_data$status == 'control'])
cor.diff.test(species_data$tspec_1[species_data$status != 'control' & species_data$longer_domain_specie == 'specie2'], species_data$tspec_2[species_data$status != 'control' & species_data$longer_domain_specie == 'specie2'],
              species_data$tspec_1[species_data$status != 'control' & species_data$longer_domain_specie == 'specie1'], species_data$tspec_2[species_data$status != 'control' & species_data$longer_domain_specie == 'specie1'])
#

#Effect of longer domain on specificity factor (ubiquitous / specificity)
c_s = sum(c(species_data$TspecF_1 == 'specific' & species_data$status != 'control' &  species_data$longer_domain_specie == 'specie2'))
c_u = sum(c(species_data$TspecF_1 == 'ubiquitous' & species_data$status != 'control' &  species_data$longer_domain_specie == 'specie2'))
m_s = sum(c(species_data$TspecF_1 == 'specific' & species_data$status != 'control' &  species_data$longer_domain_specie == 'specie1'))
m_u = sum(c(species_data$TspecF_1 == 'ubiquitous' & species_data$status != 'control' &  species_data$longer_domain_specie == 'specie1'))
 
b = barplot(c(c_s / (c_s + c_u), c_u / (c_s + c_u),
          m_s / (m_s + m_u), m_u / (m_s + m_u)),
        names.arg = rep(c('specific', 'ubiquitous'), 2),
        cex.names = 0.6, col = rep(c('blue', 'red'), each = 2), beside = T,
        main = paste0('Proportion of domain modified pair in each group of specificity factor\n', gsub('[[:digit:]]', '', species_data$GeneID_1[1])), ylab = 'pair proportion')
text(x = b, b[1,], labels = round(c(c(c_s / (c_s + c_u), c_u / (c_s + c_u)),
                           c(m_s / (m_s + m_u), m_u / (m_s + m_u))), digits = 2), cex = 0.8)
legend('bottomright', c( paste0('longest domain in ',gsub('[[:digit:]]', '', species_data$GeneID_2[1])), paste0('longest domain in ',gsub('[[:digit:]]', '', species_data$GeneID_1[1]))), fill = c('blue', 'red'), 
       cex = 0.6, horiz = F)

b = barplot(c(c_s, c_u,
              m_s, m_u),
            names.arg = rep(c('specific', 'ubiquitous'), 2),
            cex.names = 0.6, col = rep(c('blue', 'red'), each = 2), beside = T,
            main = paste0('Number of domain modified pair in each group of specificity factor\n', gsub('[[:digit:]]', '', species_data$GeneID_1[1])), ylab = 'pair proportion')
text(x = b, b[2,]*10, labels = round(c(c(c_s, c_u),
                                    c(m_s, m_u)), digits = 2), cex = 0.8)
legend('topright', c( paste0('longest domain in ',gsub('[[:digit:]]', '', species_data$GeneID_2[1])), paste0('longest domain in ',gsub('[[:digit:]]', '', species_data$GeneID_1[1]))), fill = c('blue', 'red'), 
       cex = 0.6, horiz = F)

chisq.test(species_data$TspecF_1[species_data$status != 'control'], species_data$longer_domain_specie[species_data$status != 'control'], correct = F)
####

####Paralog####
#load dataset
species_data = read.csv('D:/UNIL/Master/Master_Project/Data/expression_analysis/R_dataset/HUMAN_para_onlymale_onlybrain_dataset')
species_data$X = NULL

###Determine maximal paralog expression for a group family, according to Kryuchkova et al., 2016
#maximal expression (reference, maximal in one state) gene is the GeneID_1
#and longer domain is also GeneID_1 in modified gene
species_data = sort_paralog(paralog_dataset = species_data, expression_dataset = paste0('D:/UNIL/Master/Master_Project/Data/Bgee/', gsub('[[:digit:]]', '', species_data$GeneID_1[1]), '_expression_parsed'))
###

smoothScatter(species_data$tspec_1[species_data$status != 'control'],
              species_data$tspec_2[species_data$status != 'control'])
smoothScatter(species_data$tspec_1[species_data$status == 'control'],
              species_data$tspec_2[species_data$status == 'control'])

cor.diff.test(species_data$tspec_1[species_data$status != 'control'],
              species_data$tspec_2[species_data$status != 'control'],
              species_data$tspec_1[species_data$status == 'control'],
              species_data$tspec_2[species_data$status == 'control'])
