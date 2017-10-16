'''
Joaquim Claivaz
170919

Kaks analysis test
'''
#source("https://bioconductor.org/biocLite.R")
#biocLite("CorMut")
#biocLite("msa")
require(CorMut)
require(msa)

path_considered_ortholog = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_family/'
ortholog_family = list.files(path_considered_ortholog)
path_considered_paralog = 'D:/UNIL/Master/Master_Project/Data/OMA/paralog_family/'
paralog_family = list.files(path_considered_paralog)

#msa
mySequences = readDNAStringSet(paste0(path_considered_ortholog, ortholog_family[1]))
mySequences
myFirstAlignment = msa(mySequences)
myFirstAlignment
write.phylip(myFirstAlignment, file = 'C:/Users/Claivaz/Desktop/align')

#CorMut
# example=seqFormat(paste0(path_considered_ortholog, ortholog_family[1]), format = 'fasta')
# result=kaksCodon(example)
# fresult=filterSites(result)
# head(fresult)

s = read.alignment('C:/Users/Claivaz/Desktop/align', format = 'phylip', forceToLower = TRUE)
result = kaks(s)
as.table(result$ka/result$ks)
result$ka/result$ks
