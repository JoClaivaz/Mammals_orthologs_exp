'''
Joaquim Claivaz
171003, last modification 171012

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
####

#####Tspec analysis#####
species_data = read.table()
####Ortholog####
##DomainStatus
#Effect of domain modification (DomainStatus) on Tspec values 
hist(specie2_data$tspec[specie2_data$DomainStatus == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,specie2_data$GeneID[1])), xlab = 'Tspec value', cex.main = 0.9)
hist(specie2_data$tspec[specie2_data$DomainStatus == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')

hist(specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]), '\n', gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(specie2_data$tspec_human[specie2_data$DomainStatus == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')

#H0: no differences amongst distribution of tspec values in function of domain modification status
ks.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif'], specie2_data$tspec[specie2_data$DomainStatus == 'control'])
ks.test(specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])
#Significant: Gorgo, Macmu, Mouse, Pantr, Pigxx, Ratno
#Slighty Significant: Bovin, Mondo
#Insignificant:
#

#Effect of domain modification (DomainStatus) on specificity factor (ubiquitous / specificity)
# c_s = sum(specie2_data$TspecF[specie2_data$DomainStatus == 'control'] == 'specific')
# c_u = sum(specie2_data$TspecF[specie2_data$DomainStatus == 'control'] == 'ubiquitous')
# m_s = sum(specie2_data$TspecF[specie2_data$DomainStatus == 'modif'] == 'specific')
# m_u = sum(specie2_data$TspecF[specie2_data$DomainStatus == 'modif'] == 'ubiquitous')
# 
# barplot(c(c(c_s / (c_s + c_u), c_u / (c_s + c_u)),
#           c(m_s / (m_s + m_u), m_u / (m_s + m_u))),
#         names.arg = rep(c('specific', 'ubiquitous'), 2),
#         cex.names = 0.6, col = rep(c('blue', 'red'), each = 2), beside = T,
#         main = paste0('Pair proportion in each group of specificity factor\n', gsub('[[:digit:]]', '', specie2_data$GeneID[1])), ylab = 'pair proportion')
# legend('bottomright', c( "domain modification", "control"), fill = c('red','blue'), 
#        cex = 0.6, horiz = F)
# 
# #HO: no difference amongst specific and ubiquitous factor proportion in function of domain modification status
# chisq.test(specie2_data$TspecF, specie2_data$DomainStatus, correct = F)
#

#Correlation between tspec in function of DomainStatus
smoothScatter(specie2_data$tspec[specie2_data$DomainStatus == 'modif'], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'],
               xlab = paste0(gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' Tspec values'),
               ylab = paste0(gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]), ' Tspec values'),
               main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
                             ' and ', gsub('[[:digit:]]','' ,specie2_data$GeneID[1]),
                             '\nDomain modification group'), cex.main = 0.8)

smoothScatter(specie2_data$tspec[specie2_data$DomainStatus == 'control'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'],
               xlab = paste0(gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' Tspec values'),
               ylab = paste0(gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]), ' Tspec values'),
               main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
                             ' and ', gsub('[[:digit:]]','' ,specie2_data$GeneID[1]),
                             '\nDomain control group'), cex.main = 0.8)

# model_ancova = lm(specie2_data$tspec ~ specie2_data$tspec_human + specie2_data$DomainStatus)
# qqnorm(residuals(model_ancova)); qqline(residuals(model_ancova))
# plot(fitted.values(model_ancova), residuals(model_ancova))
# summary(model_ancova)

cor.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif'], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'])
cor.test(specie2_data$tspec[specie2_data$DomainStatus == 'control'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])

#H0: no difference in the correlation factor in tspec values pairwise comparison in function of the considered group
cor.diff.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif'], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'],
              specie2_data$tspec[specie2_data$DomainStatus == 'control'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])
#Significant: Gorgo, Mondo, Mouse, Pantr, Ratno
#Slighty Significant: Bovin
#Insignificant: Macmu, Pigxx
#

#study proportion of specific factor in function of DomainStatus
# barplot(c(table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'modif']) / sum(table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'modif'])),
#           table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'control']) / sum(table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'control']))),
#         cex.names = 0.6, las = 2, 
#         col =  rep(c('red', 'blue'), each =length(levels(specie2_data$spec_tissue))),
#         main = paste0('Pair proportion in each group of specificity factor\nIn ', gsub('[[:digit:]]', '', specie2_data$GeneID[1])), ylab = 'pair proportion')
# legend('topleft', c( "domain modification", "control"), fill = c('red','blue'), 
#        cex = 0.6, horiz = F)
# 
# barplot(c(table(specie2_data$spec_tissue_human[specie2_data$DomainStatus == 'modif']) / sum(table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'modif'])),
#           table(specie2_data$spec_tissue_human[specie2_data$DomainStatus == 'control']) / sum(table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'control']))),
#           cex.names = 0.6, col = rep(c('red', 'blue'), each =length(levels(specie2_data$spec_tissue_human))), las = 2,
#           main = paste0('Pair proportion in each group of specificity factor\nIn ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), ylab = 'pair proportion')
# legend('topleft', c( "domain modification", "control"), fill = c('red','blue'), 
#        cex = 0.6, horiz = F)
# 
# #H0: no difference amongst the proportion of specificity factor in function of domain modification
# chisq.test(specie2_data$spec_tissue, specie2_data$DomainStatus, correct = F)
# chisq.test(specie2_data$spec_tissue_human, specie2_data$DomainStatus, correct = F)
#

#study shift in specific factor in function of DomainStatus
# barplot(c(table(specie2_data$shift[specie2_data$DomainStatus == 'modif']) / sum(table(specie2_data$shift[specie2_data$DomainStatus == 'modif'])),
#           table(specie2_data$shift[specie2_data$DomainStatus == 'control']) / sum(table(specie2_data$shift[specie2_data$DomainStatus == 'control']))),
#         cex.names = 0.6, col = rep(c('red', 'blue'), each = 5), las = 2, beside = T,
#         main = paste0('Pair proportion in each group of shift specificity factor\n', gsub('[[:digit:]]', '', specie2_data$GeneID[1]), ' vs ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), ylab = 'pair proportion')
# legend('top', c( "domain modification", "control"), fill = c('red','blue'), 
#        cex = 0.6, horiz = F)

pie(table(specie2_data$shift[specie2_data$DomainStatus == 'modif']) / sum(table(specie2_data$shift[specie2_data$DomainStatus == 'modif'])),
        cex.names = 0.6 , 
        main = paste0('Pair proportion in modification group of shift specificity factor\n', gsub('[[:digit:]]', '', specie2_data$GeneID[1]), ' vs ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), ylab = 'pair proportion')

pie(table(specie2_data$shift[specie2_data$DomainStatus == 'control']) / sum(table(specie2_data$shift[specie2_data$DomainStatus == 'control'])),
    cex.names = 0.6 , 
    main = paste0('Pair proportion in control group of shift specificity factor\n', gsub('[[:digit:]]', '', specie2_data$GeneID[1]), ' vs ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), ylab = 'pair proportion')

#H0: no difference in proportion among control and modification
chisq.test(specie2_data$shift, specie2_data$DomainStatus, correct = F)

#H0: no difference in proportion of ortholog evolving in specificity ('neo- / subfunctionalization) in function of domain modification
#Shift specificity: specificity_shift, ubiquitous_to_specific
#No shift specificity: specific_no_shift, specific_to_ubiquitous, ubiquitous_no_shift
# chisq.test(c(rep(1, length(specie2_data$shift[specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift'])), 
#              rep(0, length(specie2_data$shift[!(specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift')]))), 
#            c(as.factor(specie2_data$DomainStatus[specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift']), 
#              as.factor(specie2_data$DomainStatus[!(specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift')])), correct = F)
#Significant: Mondo, Ratno
#Slighty Significant:
#Insignificant: Bovin, Gorgo, Macmu, Mouse, Pantr, Pigxx

#Shift specificity: specificity_shift, ubiquitous_to_specific, specific_to_ubiquitous
#No shift specificity: specific_no_shift, ubiquitous_no_shift
# chisq.test(c(rep(1, length(specie2_data$shift[specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift' | specie2_data$shift == 'specific_to_ubiquitous'])), 
#              rep(0, length(specie2_data$shift[!(specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift' | specie2_data$shift == 'specific_to_ubiquitous')]))), 
#            c(as.factor(specie2_data$DomainStatus[specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift' | specie2_data$shift == 'specific_to_ubiquitous']), 
#              as.factor(specie2_data$DomainStatus[!(specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift' | specie2_data$shift == 'specific_to_ubiquitous')])), correct = F)
#Significant: Mondo, Mouse, Ratno
#Slighty Significant: Pigxx
#Insignificant: Bovin, Gorgo, Macmu, Pantr
##
####

###Effect of longer domain on the different estimators####
#on Tspec values 
# hist(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], 
#      breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5), 
#      main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n', 
#                    gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' and ',gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
#                    ' comparison\nlonger domain in ', gsub('[[:digit:]]', '', specie2_data$GeneID[1])), 
#      xlab = 'Tspec value', cex.main = 0.8)
# hist(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
# legend('top', c(paste0("tspec in ", gsub('[[:digit:]]', '', specie2_data$GeneID[1])), paste0("tpsec in ", gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
# 
# hist(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], 
#      breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5), 
#      main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n', 
#                    gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' and ',gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
#                    ' comparison\nlonger domain in ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), 
#      xlab = 'Tspec value', cex.main = 0.8)
# hist(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
# legend('top', c(paste0("tspec in ", gsub('[[:digit:]]', '', specie2_data$GeneID[1])), paste0("tpsec in ", gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

# hist(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], 
#      breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5), 
#      main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n', 
#                    gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' and ',gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
#                    ' comparison\ntspec in ', gsub('[[:digit:]]', '', specie2_data$GeneID[1])), 
#      xlab = 'Tspec value', cex.main = 0.8)
# hist(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
# legend('top', c(paste0("longer domain in ", gsub('[[:digit:]]', '', specie2_data$GeneID[1])), paste0("longer domain in ", gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
# abline( v = 0.8, col = 'red')
# 
# hist(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], 
#      breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5), 
#      main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n', 
#                    gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' and ',gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
#                    ' comparison\ntspec in ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), 
#      xlab = 'Tspec value', cex.main = 0.8)
# hist(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
# legend('top', c(paste0("longer domain in ", gsub('[[:digit:]]', '', specie2_data$GeneID[1])), paste0("longer domain in ", gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
# abline( v = 0.8, col = 'red')
# 
# #H0: no difference in the tspec values distribution in function of longer domain  
# ks.test(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]) & !(is.na(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])]))],
#         specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1]) & !(is.na(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])]))])
# ks.test(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]) & !(is.na(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])]))],
#         specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1]) & !(is.na(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])]))])
# #Significant: 
# #Slighty Significant: Bovin
# #Insignificant: Gorgo, Macmu, Mondo, Mouse, Pantr, Pigxx, Ratno
# #H0: no difference in tspec value distribution between longer domain and control group
# ks.test(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])],
#         specie2_data$tspec[specie2_data$DomainStatus == 'control'])
# ks.test(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])],
#         specie2_data$tspec[specie2_data$DomainStatus == 'control'])
# ks.test(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])],
#         specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])
# ks.test(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])],
#         specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])
#Significant:
#Slighty Significant: Bovin, Gorgo, Macmu, Mouse, Pantr, Pigxx, Ratno
#Insignificant: Mondo

smoothScatter(specie2_data$tspec[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], 
              specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'  & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])],
              xlab = paste0(gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' Tspec values'),
              ylab = paste0(gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]), ' Tspec values'),
              main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
                            ' and ', gsub('[[:digit:]]','' ,specie2_data$GeneID[1]),
                            '\nDomain modification group, longer domain in ', gsub('[[:digit:]]', '', specie2_data$GeneID[1])), cex.main = 0.8)

smoothScatter(specie2_data$tspec[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], 
              specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'  & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])],
              xlab = paste0(gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' Tspec values'),
              ylab = paste0(gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]), ' Tspec values'),
              main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
                            ' and ', gsub('[[:digit:]]','' ,specie2_data$GeneID[1]),
                            '\nDomain modification group, longer domain in ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), cex.main = 0.8)

smoothScatter(specie2_data$tspec[specie2_data$DomainStatus == 'control'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'],
              xlab = paste0(gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' Tspec values'),
              ylab = paste0(gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]), ' Tspec values'),
              main = paste0('Correlation amongst Tspec values between ', gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
                            ' and ', gsub('[[:digit:]]','' ,specie2_data$GeneID[1]),
                            '\nDomain control group'), cex.main = 0.8)

cor.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'& specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])])
cor.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'& specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])])
cor.test(specie2_data$tspec[specie2_data$DomainStatus == 'control'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])

#H0: no difference in the correlation factor in tspec values pairwise comparison in function of the considered group
cor.diff.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])],
              specie2_data$tspec[specie2_data$DomainStatus == 'control'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])
cor.diff.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])],
              specie2_data$tspec[specie2_data$DomainStatus == 'control'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])
cor.diff.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])],
              specie2_data$tspec[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif' & specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])])
#Significant:
#Slighty Significant: 
#Insignificant:
#

#Effect of longer domain on specificity factor (ubiquitous / specificity)
# c_s = sum(c(specie2_data$TspecF == 'specific' & specie2_data$DomainStatus == 'modif' &  specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])))
# c_u = sum(c(specie2_data$TspecF == 'ubiquitous' & specie2_data$DomainStatus == 'modif' &  specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])))
# m_s = sum(c(specie2_data$TspecF == 'specific' & specie2_data$DomainStatus == 'modif' &  specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])))
# m_u = sum(c(specie2_data$TspecF == 'ubiquitous' & specie2_data$DomainStatus == 'modif' &  specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])))
# 
# b = barplot(c(c_s / (c_s + c_u), c_u / (c_s + c_u),
#           m_s / (m_s + m_u), m_u / (m_s + m_u)),
#         names.arg = rep(c('specific', 'ubiquitous'), 2),
#         cex.names = 0.6, col = rep(c('blue', 'red'), each = 2), beside = T,
#         main = paste0('Proportion of domain modified pair in each group of specificity factor\n', gsub('[[:digit:]]', '', specie2_data$GeneID[1])), ylab = 'pair proportion')
# text(x = b, b[1,], labels = round(c(c(c_s / (c_s + c_u), c_u / (c_s + c_u)),
#                            c(m_s / (m_s + m_u), m_u / (m_s + m_u))), digits = 2), cex = 0.8)
# legend('bottomright', c( paste0('longest domain in ',gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), paste0('longest domain in ',gsub('[[:digit:]]', '', specie2_data$GeneID[1]))), fill = c('blue', 'red'), 
#        cex = 0.6, horiz = F)
# 
# b = barplot(c(c_s, c_u,
#               m_s, m_u),
#             names.arg = rep(c('specific', 'ubiquitous'), 2),
#             cex.names = 0.6, col = rep(c('blue', 'red'), each = 2), beside = T,
#             main = paste0('Number of domain modified pair in each group of specificity factor\n', gsub('[[:digit:]]', '', specie2_data$GeneID[1])), ylab = 'pair proportion')
# text(x = b, b[2,]*10, labels = round(c(c(c_s, c_u),
#                                     c(m_s, m_u)), digits = 2), cex = 0.8)
# legend('topright', c( paste0('longest domain in ',gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), paste0('longest domain in ',gsub('[[:digit:]]', '', specie2_data$GeneID[1]))), fill = c('blue', 'red'), 
#        cex = 0.6, horiz = F)
# 
# chisq.test(specie2_data$TspecF[specie2_data$DomainStatus == 'modif'], specie2_data$longer_domain_specie[specie2_data$DomainStatus == 'modif'], correct = F)
####

####Control quality on domain modification call####
# test_data = aggregate(specie2_data$DomainStatus ~ specie2_data$GeneID + specie2_data$SpecieHomolog, FUN = paste0)
# names(test_data) = c('GeneID', 'SpecieHomolog', 'DomainStatus')
# test_data = test_data[test_data$DomainStatus == 'modif',]
# 
# path_output_pfamscan = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/'
# central_sp_name = gsub('[[:digit:]]', '', test_data$SpecieHomolog[1])
# specie2_name = gsub('[[:digit:]]', '', test_data$GeneID[1])
# 
# specie2_domain = read.table(paste0(path_output_pfamscan, specie2_name, '_domain'))
# central_sp_domain = read.table(paste0(path_output_pfamscan, central_sp_name, '_domain'))
# 
# specie2_domain = aggregate(V7 ~ V1 , data = specie2_domain, FUN = paste, collapse = '_')
# central_sp_domain = aggregate(V7 ~ V1 , data = central_sp_domain, FUN = paste, collapse = '_')
# 
# test_data$central_domain = NA
# test_data$sp2_domain = NA
# 
# for (i in 1:dim(test_data)[1]){
#   test_data$central_domain[i] = central_sp_domain$V7[central_sp_domain$V1 == as.character(test_data$SpecieHomolog[i])]
#   test_data$sp2_domain[i] = specie2_domain$V7[specie2_domain$V1 == as.character(test_data$GeneID[i])]
# }
# 
# write.csv(test_data, file = paste0('C:/Users/Claivaz/Desktop/test_domain_', specie2_name))
####

####Paralog####
smoothScatter(human_para$tspec_1[human_para$status != 'control'],
              human_para$tspec_2[human_para$status != 'control'])
smoothScatter(human_para$tspec_1[human_para$status == 'control'],
              human_para$tspec_2[human_para$status == 'control'])

cor.diff.test(human_para$tspec_1[human_para$status != 'control'],
              human_para$tspec_2[human_para$status != 'control'],
              human_para$tspec_1[human_para$status == 'control'],
              human_para$tspec_2[human_para$status == 'control'])