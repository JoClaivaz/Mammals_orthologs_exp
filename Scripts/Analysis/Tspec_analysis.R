'''
Joaquim Claivaz
171003

Tspec analysis in mammals
'''
####FUN####
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
                                     considered_sex_list = F){
  
  specie2_data = read.table(paste0('D:/UNIL/Master/Master_Project/Data/Bgee/',
                                   considered_species_name,
                                   '_expression_parsed'), header = T, sep = '\t')
  central_sp_data = read.table(paste0('D:/UNIL/Master/Master_Project/Data/Bgee/',
                                      central_species_name,
                                      '_expression_parsed'), header = T, sep = '\t')
  
  #Sex
  if (considered_sex_list != F){
    central_sp_data = central_sp_data[central_sp_data$Sex %in% considered_sex_list,]
    specie2_data = specie2_data[specie2_data$Sex %in% considered_sex_list,]
    
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

data_organization_sex_consideration = function(considered_species_name, central_species_name = 'HUMAN'){
  
  #DOESN'T WORK WELL, NEED CONSIDERATION OF SOME NA IN SEX FOR SOME SPECIES (such as replace NA by 'undefined')
  
  specie2_data = read.table(paste0('D:/UNIL/Master/Master_Project/Data/Bgee/',
                                   considered_species_name,
                                   '_expression_parsed'), header = T, sep = '\t')
  central_sp_data = read.table(paste0('D:/UNIL/Master/Master_Project/Data/Bgee/',
                                      central_species_name,
                                      '_expression_parsed'), header = T, sep = '\t')
  
  #Sex
  notkeep_Sex = unique(central_sp_data$Sex)[!(unique(central_sp_data$Sex) %in% unique(specie2_data$Sex))]
  central_sp_data = central_sp_data[!(central_sp_data$Sex %in% notkeep_Sex),]
  specie2_data = specie2_data[!(specie2_data$Sex %in% notkeep_Sex),]
  #
  
  ##intersect of dataset
  #AnatomicalEntityName
  keep_AnatomicalEntityName = unique(specie2_data$AnatomicalEntityName)[unique(specie2_data$AnatomicalEntityName) %in% unique(central_sp_data$AnatomicalEntityName)]
  central_sp_data = central_sp_data[central_sp_data$AnatomicalEntityName %in% keep_AnatomicalEntityName,]
  specie2_data = specie2_data[specie2_data$AnatomicalEntityName %in% keep_AnatomicalEntityName,]
  #
  
  #MEAN fpkm in function of geneID and AnatomicalEntityName
  specie2_data = aggregate(specie2_data$FPKM ~ specie2_data$GeneID + specie2_data$AnatomicalEntityName + specie2_data$HUMAN_homolog + specie2_data$DomainStatus + specie2_data$Sex, FUN = mean)
  names(specie2_data) = c('GeneID', 'Tissue', 'SpecieHomolog', 'DomainStatus', 'Sex', 'FPKM')
  central_sp_data = aggregate(central_sp_data$FPKM ~ central_sp_data$GeneID + central_sp_data$AnatomicalEntityName + central_sp_data$Sex, FUN = mean)
  names(central_sp_data) = c('GeneID', 'Tissue', 'Sex', 'FPKM')
  ##
  
  #Group data
  all(unique(as.character(specie2_data$Tissue)) == unique(as.character(central_sp_data$Tissue)))
  names_tissue = unique(as.character(specie2_data$Tissue))
  central_sp_data = aggregate(FPKM ~ GeneID + Sex, data = central_sp_data, paste, collapse = '_')
  central_sp_data = separate(central_sp_data, col = 'FPKM', sep= '_', into = names_tissue)
  specie2_data = aggregate(FPKM ~ GeneID + SpecieHomolog + DomainStatus + Sex, data = specie2_data, paste, collapse = '_')
  specie2_data = separate(specie2_data, col = 'FPKM', sep= '_', into = names_tissue)
  #
  
  #convert fpkm into numerics
  for (col_tmp in 3:dim(central_sp_data)[2]){
    central_sp_data[,col_tmp] = as.numeric(central_sp_data[,col_tmp])
  }
  
  for (col_tmp in 5:dim(specie2_data)[2]){
    specie2_data[,col_tmp] = as.numeric(specie2_data[,col_tmp])
  }
  #
  
  #Tspec inference with log_transformation
  central_sp_data = Tspec_inference(central_sp_data, first_column_num = 3)
  specie2_data = Tspec_inference(specie2_data, first_column_num = 5)
  #
  
  #specificity factor inference
  central_sp_data = specificity_inference(central_sp_data, first_column_num = 3)
  specie2_data = specificity_inference(specie2_data, first_column_num = 5)
  #
  
  #dataframe containing Specie2 and central species data
  col_keep = c('GeneID', 'SpecieHomolog', 'DomainStatus', 'Sex', 'tspec', 'TspecF', 'spec_tissue')
  specie2_data = specie2_data[,names(specie2_data) %in% col_keep]
  central_sp_data = central_sp_data[,names(central_sp_data) %in% col_keep]
  names(central_sp_data) = c('GeneID_human', 'Sex', 'tspec_human', 'TspecF_human', 'spec_tissue_human') 
  
  specie2_data$tspec_human = NA
  specie2_data$TspecF_human = NA
  specie2_data$spec_tissue_human = NA
  
  for (gene_row in 1:dim(specie2_data)[1]){
    row_central_sp = which(as.character(central_sp_data$GeneID_human) == as.character(specie2_data$SpecieHomolog[gene_row])
                           )[which(as.character(central_sp_data$GeneID_human) == as.character(specie2_data$SpecieHomolog[gene_row])
                           ) %in% which(as.character(central_sp_data$Sex) == as.character(specie2_data$Sex[gene_row]))]
    
    if(length(central_sp_data$tspec_human[row_central_sp]) > 0){
      specie2_data$tspec_human[gene_row] = central_sp_data$tspec_human[row_central_sp]
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
  specie2_data$Sex = as.factor(specie2_data$Sex)
  specie2_data$tspec = as.numeric(specie2_data$tspec)
  specie2_data$tspec_human = as.numeric(specie2_data$tspec_human)
  
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
                                       considered_sex_list = c('male'))

#
# gorgo_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/GORGO_expression_parsed', header = T, sep = '\t')
# unique(gorgo_test$AnatomicalEntityName)[unique(gorgo_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(gorgo_test$StageName)[unique(gorgo_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(gorgo_test$Sex)[unique(gorgo_test$Sex) %in% unique(human_test$Sex)]
#need sex consideration (female / male)
gorgo_data = data_organization_one_sex(considered_species_name = 'GORGO',
                                       considered_sex_list = c('male'))

#
# macmu_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/MACMU_expression_parsed', header = T, sep = '\t')
# unique(macmu_test$AnatomicalEntityName)[unique(macmu_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(macmu_test$StageName)[unique(macmu_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(macmu_test$Sex)[unique(macmu_test$Sex) %in% unique(human_test$Sex)]
#sex consideration (female / male / NA)
macmu_data = data_organization_one_sex(considered_species_name = 'MACMU',
                                       considered_sex_list = c('male'))

#
# mondo_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/MONDO_expression_parsed', header = T, sep = '\t')
# unique(mondo_test$AnatomicalEntityName)[unique(mondo_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(mondo_test$StageName)[unique(mondo_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(mondo_test$Sex)[unique(mondo_test$Sex) %in% unique(human_test$Sex)]
#sex consideration (female / male)
mondo_data = data_organization_one_sex(considered_species_name = 'MONDO',
                                       considered_sex_list = c('male'))

#
# mouse_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/MOUSE_expression_parsed', header = T, sep = '\t')
# unique(mouse_test$AnatomicalEntityName)[unique(mouse_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(mouse_test$StageName)[unique(mouse_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(mouse_test$Sex)[unique(mouse_test$Sex) %in% unique(human_test$Sex)]
#sex consideration !warning sex factor: female male NA mixed (No consider mixed, as first step)
mouse_data = data_organization_one_sex(considered_species_name = 'MOUSE',
                                       considered_sex_list = c('male'))

#
# pantr_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/PANTR_expression_parsed', header = T, sep = '\t')
# unique(pantr_test$AnatomicalEntityName)[unique(pantr_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(pantr_test$StageName)[unique(pantr_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(pantr_test$Sex)[unique(pantr_test$Sex) %in% unique(human_test$Sex)]
#sex consideration (female / male)
pantr_data = data_organization_one_sex(considered_species_name = 'PANTR',
                                       considered_sex_list = c('male'))

#
# pigxx_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/pigxx_expression_parsed', header = T, sep = '\t')
# unique(pigxx_test$AnatomicalEntityName)[unique(pigxx_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(pigxx_test$StageName)[unique(pigxx_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(pigxx_test$Sex)[unique(pigxx_test$Sex) %in% unique(human_test$Sex)]
#as BOVIN
pigxx_data = data_organization_one_sex(considered_species_name = 'PIGXX',
                                       considered_sex_list = c('male'))

#
# ratno_test = read.table('D:/UNIL/Master/Master_Project/Data/Bgee/RATNO_expression_parsed', header = T, sep = '\t')
# unique(ratno_test$AnatomicalEntityName)[unique(ratno_test$AnatomicalEntityName)
#                                         %in% unique(human_test$AnatomicalEntityName)]
# unique(ratno_test$StageName)[unique(ratno_test$StageName)
#                              %in% unique(human_test$StageName)]
# unique(ratno_test$Sex)[unique(ratno_test$Sex) %in% unique(human_test$Sex)]
#as BOVIN
ratno_data = data_organization_one_sex(considered_species_name = 'RATNO',
                                       considered_sex_list = c('male'))
####

####Tspec analysis####
##Dataset available / consider only male analysis
specie2_data = bovin_data
specie2_data = gorgo_data
specie2_data = macmu_data
specie2_data = mondo_data
specie2_data = mouse_data
specie2_data = pantr_data
specie2_data = pigxx_data
specie2_data = ratno_data
#Distribution effect of domain modification (DomainStatus) on Tspec values 
hist(specie2_data$tspec[specie2_data$DomainStatus == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,specie2_data$GeneID[1])), xlab = 'Tspec value', cex.main = 0.9)
hist(specie2_data$tspec[specie2_data$DomainStatus == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(specie2_data$tspec[specie2_data$DomainStatus == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), 
     main = paste0('Distribution of tissue specificity values in ', gsub('[[:digit:]]','' ,specie2_data$SpecieHomolog[1]), '\n', gsub('[[:digit:]]','' ,specie2_data$GeneID[1]), ' comparison'), xlab = 'Tspec value', cex.main = 0.9)
hist(specie2_data$tspec[specie2_data$DomainStatus == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

kruskal.test(specie2_data$tspec ~ specie2_data$DomainStatus)
kruskal.test(specie2_data$tspec_human ~ specie2_data$DomainStatus)
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

model_ancova = lm(specie2_data$tspec ~ specie2_data$tspec_human + specie2_data$DomainStatus)
qqnorm(residuals(model_ancova)); qqline(residuals(model_ancova))
plot(fitted.values(model_ancova), residuals(model_ancova))
summary(model_ancova)

cor.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif'], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'])
cor.test(specie2_data$tspec[specie2_data$DomainStatus == 'control'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])

cor.diff.test(specie2_data$tspec[specie2_data$DomainStatus == 'modif'], specie2_data$tspec_human[specie2_data$DomainStatus == 'modif'],
              specie2_data$tspec[specie2_data$DomainStatus == 'control'], specie2_data$tspec_human[specie2_data$DomainStatus == 'control'])
#better results for BOVIN, if male and NA are used using: 
#data_organization_one_sex(considered_species_name = 'BOVIN', considered_sex_list = F)

#study proportion of specific factor in function of DomainStatus
barplot(c(table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'modif']) / sum(table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'modif'])),
          table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'control']) / table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'control'])),
        cex.names = 0.6, col = c('red', 'blue'), las = 2,
        names.arg = rep(names(table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'control'])), each = 2),
        main = paste0('Pair proportion in each group of specificity factor\n', gsub('[[:digit:]]', '', specie2_data$GeneID[1])), ylab = 'pair proportion')
legend('top', c( "domain modification", "control"), fill = c('red','blue'), 
       cex = 0.6, horiz = F)

barplot(c(table(specie2_data$spec_tissue_human[specie2_data$DomainStatus == 'modif']) / sum(table(specie2_data$spec_tissue_human[specie2_data$DomainStatus == 'modif'])),
          table(specie2_data$spec_tissue_human[specie2_data$DomainStatus == 'control']) / sum(table(specie2_data$spec_tissue_human[specie2_data$DomainStatus == 'control']))),
        cex.names = 0.6, col = c('red', 'blue'), las = 2,
        names.arg = rep(names(table(specie2_data$spec_tissue[specie2_data$DomainStatus == 'control'])), each = 2),
        main = paste0('Pair proporton in each group of specificity factor\n', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), ylab = 'pair proportion')
legend('top', c( "domain modification", "control"), fill = c('red','blue'), 
       cex = 0.6, horiz = F)

chisq.test(specie2_data$spec_tissue, specie2_data$DomainStatus, correct = F)

chisq.test(specie2_data$spec_tissue_human, specie2_data$DomainStatus, correct = F)
#

#study shift in specific factor in function of DomainStatus
barplot(c(table(specie2_data$shift[specie2_data$DomainStatus == 'modif']) / sum(table(specie2_data$shift[specie2_data$DomainStatus == 'modif'])),
          table(specie2_data$shift[specie2_data$DomainStatus == 'control']) / table(specie2_data$shift[specie2_data$DomainStatus == 'control'])),
        cex.names = 0.6, col = c('red', 'blue'), las = 2,
        names.arg = rep(names(table(specie2_data$shift[specie2_data$DomainStatus == 'control'])), each = 2),
        main = paste0('Proportion pair in each group of shift specificity factor\n', gsub('[[:digit:]]', '', specie2_data$GeneID[1]), ' vs ', gsub('[[:digit:]]', '', specie2_data$SpecieHomolog[1])), ylab = 'pair proportion')
legend('top', c( "domain modification", "control"), fill = c('red','blue'), 
       cex = 0.6, horiz = F)

chisq.test(specie2_data$shift ,specie2_data$DomainStatus, correct = F)
##
####