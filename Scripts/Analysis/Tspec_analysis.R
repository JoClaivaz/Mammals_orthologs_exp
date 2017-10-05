'''
Joaquim Claivaz
171003

Tspec analysis in mammals
'''
####FUN####
##Tspec analysis##
log_transformation_rpkm = function(data_frame, first_column_num){
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

Tspec_inference = function(expression_data, first_column_num){
  expression_data = log_transformation_rpkm(expression_data, first_column_num)
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

data_organization_one_sex = function(considered_species_name, 
                                     central_species_name = 'HUMAN',
                                     considered_sex_vector = F,
                                     considered_anat_vector = F,
                                     considered_devtime_vector = F){
  
  specie2_data = read.table(paste0('D:/UNIL/Master/Master_Project/Data/Bgee/',
                                   considered_species_name,
                                   '_expression_parsed'), header = T, sep = '\t')
  central_sp_data = read.table(paste0('D:/UNIL/Master/Master_Project/Data/Bgee/',
                                      central_species_name,
                                      '_expression_parsed'), header = T, sep = '\t')
  
  #Sex
  if (considered_sex_vector != F){
    central_sp_data = central_sp_data[central_sp_data$Sex %in% considered_sex_vector,]
    specie2_data = specie2_data[specie2_data$Sex %in% considered_sex_vector,]
    
  }else{
    notkeep_Sex = unique(central_sp_data$Sex)[!(unique(central_sp_data$Sex) %in% unique(specie2_data$Sex))]
    central_sp_data = central_sp_data[!(central_sp_data$Sex %in% notkeep_Sex),]
    specie2_data = specie2_data[!(specie2_data$Sex %in% notkeep_Sex),]
  }
  #
  
  ##intersect of dataset
  #AnatomicalEntityName
  keep_AnatomicalEntityName = unique(specie2_data$AnatomicalEntityName)[unique(specie2_data$AnatomicalEntityName) %in% unique(central_sp_data$AnatomicalEntityName)]
  central_sp_data = central_sp_data[central_sp_data$AnatomicalEntityName %in% keep_AnatomicalEntityName,]
  specie2_data = specie2_data[specie2_data$AnatomicalEntityName %in% keep_AnatomicalEntityName,]
  #
  
  if (considered_anat_vector != F){
    central_sp_data = central_sp_data[central_sp_data$AnatomicalEntityName %in% considered_anat_vector,]
    specie2_data = specie2_data[specie2_data$AnatomicalEntityName %in% considered_anat_vector,]
    
  }
  
  #developmental time filtering
  if (considered_devtime_vector != F){
    central_sp_data = central_sp_data[central_sp_data$StageName %in% considered_devtime_vector,]
    specie2_data = specie2_data[specie2_data$StageName %in% considered_devtime_vector,]
    
  }
  #
  
  #MEAN fpkm in function of geneID and AnatomicalEntityName
  specie2_data = aggregate(specie2_data$FPKM ~ specie2_data$GeneID + specie2_data$AnatomicalEntityName + specie2_data$HUMAN_homolog + specie2_data$DomainStatus, FUN = mean)
  names(specie2_data) = c('GeneID', 'Tissue', 'SpecieHomolog', 'DomainStatus', 'FPKM')
  central_sp_data = aggregate(central_sp_data$FPKM ~ central_sp_data$GeneID + central_sp_data$AnatomicalEntityName, FUN = mean)
  names(central_sp_data) = c('GeneID', 'Tissue', 'FPKM')
  ##
  
  #Group data
  all(unique(as.character(specie2_data$Tissue)) == unique(as.character(central_sp_data$Tissue)))
  names_tissue = unique(as.character(specie2_data$Tissue))
  central_sp_data = aggregate(FPKM ~ GeneID, data = central_sp_data, paste, collapse = '_')
  central_sp_data = separate(central_sp_data, col = 'FPKM', sep= '_', into = names_tissue)
  specie2_data = aggregate(FPKM ~ GeneID + SpecieHomolog + DomainStatus, data = specie2_data, paste, collapse = '_')
  specie2_data = separate(specie2_data, col = 'FPKM', sep= '_', into = names_tissue)
  #
  
  #convert fpkm into numerics
  for (col_tmp in 2:dim(central_sp_data)[2]){
    central_sp_data[,col_tmp] = as.numeric(central_sp_data[,col_tmp])
  }
  
  for (col_tmp in 4:dim(specie2_data)[2]){
    specie2_data[,col_tmp] = as.numeric(specie2_data[,col_tmp])
  }
  #
  
  #Tspec inference with log_transformation
  central_sp_data = Tspec_inference(central_sp_data, first_column_num = 2)
  specie2_data = Tspec_inference(specie2_data, first_column_num = 4)
  #
  
  #specificity factor inference
  central_sp_data = specificity_inference(central_sp_data, first_column_num = 2)
  specie2_data = specificity_inference(specie2_data, first_column_num = 4)
  #
  
  #dataframe containing Specie2 and central species data
  col_keep = c('GeneID', 'SpecieHomolog', 'DomainStatus', 'tspec', 'TspecF', 'spec_tissue')
  specie2_data = specie2_data[,names(specie2_data) %in% col_keep]
  central_sp_data = central_sp_data[,names(central_sp_data) %in% col_keep]
  names(central_sp_data) = c('GeneID_human', 'tspec_human', 'TspecF_human', 'spec_tissue_human') 
  
  specie2_data$tspec_human = NA
  specie2_data$TspecF_human = NA
  specie2_data$spec_tissue_human = NA
  
  for (gene_row in 1:dim(specie2_data)[1]){
    row_central_sp = which(as.character(central_sp_data$GeneID_human) == as.character(specie2_data$SpecieHomolog[gene_row]))
    
    if(length(central_sp_data$tspec_human[row_central_sp]) > 0){
      specie2_data$tspec_human[gene_row] = as.character(central_sp_data$tspec_human[row_central_sp])
      specie2_data$TspecF_human[gene_row] = as.character(central_sp_data$TspecF_human[row_central_sp])
      specie2_data$spec_tissue_human[gene_row] = as.character(central_sp_data$spec_tissue_human[row_central_sp])
    }
  }
  #
  
  #keep only complete case
  specie2_data = specie2_data[complete.cases(specie2_data),]
  #
  
  #specificity shift inference
  specie2_data$shift = NA
  for (pair_row in 1:dim(specie2_data)[1]){
    if (as.character(specie2_data$spec_tissue[pair_row]) == as.character(specie2_data$spec_tissue_human[pair_row])){
      if(as.character(specie2_data$spec_tissue[pair_row]) != 'ubiquitous'){
        specie2_data$shift[pair_row] = 'specific_no_shift'
      }else{
        specie2_data$shift[pair_row] = 'ubiquitous_no_shift'
      }
    }
    else if (all(as.character(specie2_data$spec_tissue[pair_row]) != 'ubiquitous',
                 as.character(specie2_data$spec_tissue_human[pair_row]) != 'ubiquitous')){
      specie2_data$shift[pair_row] = 'specificity_shift'
    }
    else if (as.character(specie2_data$spec_tissue[pair_row]) != 'ubiquitous'){
      specie2_data$shift[pair_row] = 'specific_to_ubiquitous'
    }
    else if (as.character(specie2_data$spec_tissue_human[pair_row]) != 'ubiquitous'){
      specie2_data$shift[pair_row] = 'ubiquitous_to_specific'
    }
    else{
      specie2_data$shift[pair_row] = 'undefined'
    }
  }
  
  specie2_data$DomainStatus = as.factor(specie2_data$DomainStatus)
  specie2_data$TspecF = as.factor(specie2_data$TspecF)
  specie2_data$spec_tissue = as.factor(specie2_data$spec_tissue)
  specie2_data$TspecF_human = as.factor(specie2_data$TspecF_human)
  specie2_data$spec_tissue_human = as.factor(specie2_data$spec_tissue_human)
  specie2_data$shift = as.factor(specie2_data$shift)
  specie2_data$tspec = as.numeric(specie2_data$tspec)
  specie2_data$tspec_human = as.numeric(specie2_data$tspec_human)
  
  return(specie2_data)
}
##

##Effect of longer domain
##
add_specie_longer_domain_factor = function(specie2_data, 
                                           path_output_pfamscan = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/'){
  central_sp_name = gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])
  specie2_name = gsub('[[:digit:]]', '', specie2_data$GeneID[1])
  
  specie2_domain = read.table(paste0(path_output_pfamscan, specie2_name, '_domain'))
  central_sp_domain = read.table(paste0(path_output_pfamscan, central_sp_name, '_domain'))
  
  specie2_domain = aggregate(V7 ~ V1 , data = specie2_domain, FUN = length)
  central_sp_domain = aggregate(V7 ~ V1 , data = central_sp_domain, FUN = length)
  
  specie2_data$length_specie2 = NA
  specie2_data$length_central_sp = NA
  specie2_data$longer_domain_specie = NA
  
  for (gene_row in 1:dim(specie2_data)[1]){
    if (specie2_data$DomainStatus[gene_row] == 'modif'){
      specie2_data$length_specie2[gene_row] = specie2_domain$V7[specie2_domain$V1 == as.character(specie2_data$GeneID[gene_row])]
      specie2_data$length_central_sp[gene_row] = central_sp_domain$V7[central_sp_domain$V1 == as.character(specie2_data$SpecieHomolog[gene_row])]
      
      if (specie2_data$length_specie2[gene_row] > specie2_data$length_central_sp[gene_row]){
        specie2_data$longer_domain_specie[gene_row] = specie2_name
      }
      else{
        specie2_data$longer_domain_specie[gene_row] = central_sp_name 
      }
    }
  }
  
  specie2_data$longer_domain_specie = as.factor(specie2_data$longer_domain_specie)
  
  return(specie2_data)
}
####

####Library needed####
require(tidyr)
####

####Data organization####
#human_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/HUMAN_expression_parsed', header = T, sep = '\t')
#
# bovin_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/BOVIN_expression_parsed', header = T, sep = '\t')
# unique(bovin_test$AnatomicalEntityName)[unique(bovin_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(bovin_test$StageName)[unique(bovin_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(bovin_test$Sex)[unique(bovin_test$Sex) %in% unique(human_test$Sex)]
#sex available (male / NA)
bovin_data = data_organization_one_sex(considered_species_name = 'BOVIN',
                                       considered_sex_vector = c('male'))

#
# gorgo_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/GORGO_expression_parsed', header = T, sep = '\t')
# unique(gorgo_test$AnatomicalEntityName)[unique(gorgo_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(gorgo_test$StageName)[unique(gorgo_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(gorgo_test$Sex)[unique(gorgo_test$Sex) %in% unique(human_test$Sex)]
#need sex consideration (female / male)
gorgo_data = data_organization_one_sex(considered_species_name = 'GORGO',
                                       considered_sex_vector = c('male'))

#
# macmu_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/MACMU_expression_parsed', header = T, sep = '\t')
# unique(macmu_test$AnatomicalEntityName)[unique(macmu_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(macmu_test$StageName)[unique(macmu_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(macmu_test$Sex)[unique(macmu_test$Sex) %in% unique(human_test$Sex)]
#sex consideration (female / male / NA)
macmu_data = data_organization_one_sex(considered_species_name = 'MACMU',
                                       considered_sex_vector = c('male'))

#
# mondo_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/MONDO_expression_parsed', header = T, sep = '\t')
# unique(mondo_test$AnatomicalEntityName)[unique(mondo_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(mondo_test$StageName)[unique(mondo_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(mondo_test$Sex)[unique(mondo_test$Sex) %in% unique(human_test$Sex)]
#sex consideration (female / male)
mondo_data = data_organization_one_sex(considered_species_name = 'MONDO',
                                       considered_sex_vector = c('male'))

#
# mouse_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/MOUSE_expression_parsed', header = T, sep = '\t')
# unique(mouse_test$AnatomicalEntityName)[unique(mouse_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(mouse_test$StageName)[unique(mouse_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(mouse_test$Sex)[unique(mouse_test$Sex) %in% unique(human_test$Sex)]
#sex consideration !warning sex factor: female male NA mixed (No consider mixed, as first step)
mouse_data = data_organization_one_sex(considered_species_name = 'MOUSE',
                                       considered_sex_vector = c('male'))

#
# pantr_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/PANTR_expression_parsed', header = T, sep = '\t')
# unique(pantr_test$AnatomicalEntityName)[unique(pantr_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(pantr_test$StageName)[unique(pantr_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(pantr_test$Sex)[unique(pantr_test$Sex) %in% unique(human_test$Sex)]
#sex consideration (female / male)
pantr_data = data_organization_one_sex(considered_species_name = 'PANTR',
                                       considered_sex_vector = c('male'))

#
# pigxx_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/pigxx_expression_parsed', header = T, sep = '\t')
# unique(pigxx_test$AnatomicalEntityName)[unique(pigxx_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(pigxx_test$StageName)[unique(pigxx_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(pigxx_test$Sex)[unique(pigxx_test$Sex) %in% unique(human_test$Sex)]
#as BOVIN
pigxx_data = data_organization_one_sex(considered_species_name = 'PIGXX',
                                       considered_sex_vector = c('male'))

#
# ratno_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/RATNO_expression_parsed', header = T, sep = '\t')
# unique(ratno_test$AnatomicalEntityName)[unique(ratno_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(ratno_test$StageName)[unique(ratno_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(ratno_test$Sex)[unique(ratno_test$Sex) %in% unique(human_test$Sex)]
#as BOVIN
ratno_data = data_organization_one_sex(considered_species_name = 'RATNO',
                                       considered_sex_vector = c('male'))

#add domain architecture length and infer which species has the bigger domain
bovin_data = add_specie_longer_domain_factor(bovin_data)
gorgo_data = add_specie_longer_domain_factor(gorgo_data)
macmu_data = add_specie_longer_domain_factor(macmu_data)
mondo_data = add_specie_longer_domain_factor(mondo_data)
mouse_data = add_specie_longer_domain_factor(mouse_data)
pantr_data = add_specie_longer_domain_factor(pantr_data)
pigxx_data = add_specie_longer_domain_factor(pigxx_data)
ratno_data = add_specie_longer_domain_factor(ratno_data)
####

####Tspec analysis####
###Dataset available / considered only male analysis
specie2_data = bovin_data
specie2_data = gorgo_data
specie2_data = macmu_data
specie2_data = mondo_data
specie2_data = mouse_data
specie2_data = pantr_data
specie2_data = pigxx_data
specie2_data = ratno_data

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
barplot(c(table(specie2_data$shift[specie2_data$DomainStatus == 'modif']) / sum(table(specie2_data$shift[specie2_data$DomainStatus == 'modif'])),
          table(specie2_data$shift[specie2_data$DomainStatus == 'control']) / sum(table(specie2_data$shift[specie2_data$DomainStatus == 'control']))),
        cex.names = 0.6, col = rep(c('red', 'blue'), each = 5), las = 2, beside = T,
        main = paste0('Pair proportion in each group of shift specificity factor\n', gsub('[[:digit:]]', '', specie2_data$GeneID[1]), ' vs ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), ylab = 'pair proportion')
legend('top', c( "domain modification", "control"), fill = c('red','blue'), 
       cex = 0.6, horiz = F)

#H0: no difference in proportion of ortholog evolving in specificity ('neo- / subfunctionalization) in function of domain modification
#Shift specificity: specificity_shift, ubiquitous_to_specific
#No shift specificity: specific_no_shift, specific_to_ubiquitous, ubiquitous_no_shift
chisq.test(c(rep(1, length(specie2_data$shift[specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift'])), 
             rep(0, length(specie2_data$shift[!(specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift')]))), 
           c(as.factor(specie2_data$DomainStatus[specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift']), 
             as.factor(specie2_data$DomainStatus[!(specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift')])), correct = F)
#Significant: Mondo, Ratno
#Slighty Significant:
#Insignificant: Bovin, Gorgo, Macmu, Mouse, Pantr, Pigxx

#Shift specificity: specificity_shift, ubiquitous_to_specific, specific_to_ubiquitous
#No shift specificity: specific_no_shift, ubiquitous_no_shift
chisq.test(c(rep(1, length(specie2_data$shift[specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift' | specie2_data$shift == 'specific_to_ubiquitous'])), 
             rep(0, length(specie2_data$shift[!(specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift' | specie2_data$shift == 'specific_to_ubiquitous')]))), 
           c(as.factor(specie2_data$DomainStatus[specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift' | specie2_data$shift == 'specific_to_ubiquitous']), 
             as.factor(specie2_data$DomainStatus[!(specie2_data$shift == 'ubiquitous_to_specific' | specie2_data$shift == 'specificity_shift' | specie2_data$shift == 'specific_to_ubiquitous')])), correct = F)
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

hist(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], 
     breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n', 
                   gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' and ',gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
                   ' comparison\ntspec in ', gsub('[[:digit:]]', '', specie2_data$GeneID[1])), 
     xlab = 'Tspec value', cex.main = 0.8)
hist(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c(paste0("longer domain in ", gsub('[[:digit:]]', '', specie2_data$GeneID[1])), paste0("longer domain in ", gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')

hist(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])], 
     breaks = 20, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Domain modification group tissue specificity values distribution in function of domain length\n', 
                   gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' and ',gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]),
                   ' comparison\ntspec in ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), 
     xlab = 'Tspec value', cex.main = 0.8)
hist(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])], breaks = 20, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c(paste0("longer domain in ", gsub('[[:digit:]]', '', specie2_data$GeneID[1])), paste0("longer domain in ", gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]))), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
abline( v = 0.8, col = 'red')

#H0: no difference in the tspec values distribution in function of longer domain  
ks.test(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]) & !(is.na(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])]))],
        specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1]) & !(is.na(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])]))])
ks.test(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1]) & !(is.na(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])]))],
        specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1]) & !(is.na(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])]))])
#Significant: 
#Slighty Significant: Bovin
#Insignificant: Gorgo, Macmu, Mondo, Mouse, Pantr, Pigxx, Ratno
#H0: no difference in tspec value distribution between longer domain and control group
ks.test(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])],
        specie2_data$tspec[specie2_data$DomainStatus == 'control'])
ks.test(specie2_data$tspec[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])],
        specie2_data$tspec[specie2_data$DomainStatus == 'control'])
ks.test(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])],
        specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])
ks.test(specie2_data$tspec_human[specie2_data$longer_domain_specie == gsub('[[:digit:]]', '', specie2_data$GeneID[1])],
        specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])
#Significant:
#Slighty Significant: Bovin, Gorgo, Macmu, Mouse, Pantr, Pigxx, Ratno
#Insignificant: Mondo
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