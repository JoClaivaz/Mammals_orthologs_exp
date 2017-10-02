'''
Joaquim Claivaz
171001

plot rsults of domain rearrangement in mammals analysis
'''

""" dict output from analysis python
{'BOVIN': {'OMA': 15591,
  'domainDIFF': 224786,
'final': 352,
'group': {'1_domain_modif': 1056, 'complex_modif': 462, 'nomodif': 13087},
'pfamscan': {'both': 14605, 'none': 822, 'partial': 164},
'repeat': {'1_domain_complex_modif': 50,
'1_domain_notrepeated': 472,
'1_domain_repeated': 534}},
'GORGO': {'OMA': 16349,
'domainDIFF': 232346,
'final': 319,
'group': {'1_domain_modif': 896, 'complex_modif': 323, 'nomodif': 13979},
'pfamscan': {'both': 15198, 'none': 986, 'partial': 165},
'repeat': {'1_domain_complex_modif': 16,
'1_domain_notrepeated': 399,
'1_domain_repeated': 481}},
'MACMU': {'OMA': 15438,
'domainDIFF': 257597,
'final': 374,
'group': {'1_domain_modif': 1024, 'complex_modif': 351, 'nomodif': 13043},
'pfamscan': {'both': 14418, 'none': 830, 'partial': 190},
'repeat': {'1_domain_complex_modif': 27,
'1_domain_notrepeated': 480,
'1_domain_repeated': 517}},
'MONDO': {'OMA': 10564,
'domainDIFF': 183774,
'final': 319,
'group': {'1_domain_modif': 926, 'complex_modif': 426, 'nomodif': 8555},
'pfamscan': {'both': 9907, 'none': 492, 'partial': 165},
'repeat': {'1_domain_complex_modif': 42,
'1_domain_notrepeated': 424,
'1_domain_repeated': 460}},
'MOUSE': {'OMA': 16051,
'domainDIFF': 272425,
'final': 276,
'group': {'1_domain_modif': 923, 'complex_modif': 391, 'nomodif': 13676},
'pfamscan': {'both': 14990, 'none': 892, 'partial': 169},
'repeat': {'1_domain_complex_modif': 42,
'1_domain_notrepeated': 392,
'1_domain_repeated': 489}},
'PANTR': {'OMA': 16447,
'domainDIFF': 206993,
'final': 242,
'group': {'1_domain_modif': 734, 'complex_modif': 246, 'nomodif': 14353},
'pfamscan': {'both': 15333, 'none': 994, 'partial': 120},
'repeat': {'1_domain_complex_modif': 12,
'1_domain_notrepeated': 300,
'1_domain_repeated': 422}},
'PIGXX': {'OMA': 12485,
'domainDIFF': 228663,
'final': 426,
'group': {'1_domain_modif': 1010, 'complex_modif': 424, 'nomodif': 10214},
'pfamscan': {'both': 11648, 'none': 647, 'partial': 190},
'repeat': {'1_domain_complex_modif': 36,
'1_domain_notrepeated': 544,
'1_domain_repeated': 430}},
'RATNO': {'OMA': 15490,
'domainDIFF': 242116,
'final': 354,
'group': {'1_domain_modif': 1041, 'complex_modif': 446, 'nomodif': 12993},
'pfamscan': {'both': 14480, 'none': 824, 'partial': 186},
'repeat': {'1_domain_complex_modif': 42,
'1_domain_notrepeated': 473,
'1_domain_repeated': 526}}}

N.B.: complex_modif and 1_domain_complex_modif are taken as the same group
Red group represents always the group of interest
"""

#number of unique ortholog called from OMA
dat = read.table(text = "BOVIN GORGO MACMU MONDO MOUSE PANTR PIGXX RATNO
1 15591 16349 15438 10564 16051 16447 12485 15490", header = T)
b = barplot(as.matrix(dat), col = c('red'), ylim = c(0, 17000),
            ylab = 'Pair count', main = 'Number of unique ortholog called from OMA',
            cex.main = 0.7, cex.names = 0.8)
text(x = b, dat[1,] / 2, labels = dat[], cex = 0.8)

#pfamscan domain inferred by ortholog pair (both/partial/none)
dat = read.table(text = "BOVIN GORGO MACMU MONDO MOUSE PANTR PIGXX RATNO
1 14605 15198 14418 9907 14990 15333 11648 14484
2 822   986   830   492  892   994   647   824
3 164   165   190   165  169   120   190   186", header = T)
b = barplot(as.matrix(dat), col = c('red', 'blue', 'green'),
            ylab = 'Pair count', main = 'pfamscan domain inferred by ortholog pair',
            cex.main = 0.7, cex.names = 0.8, ylim = c(0, max(dat) + 3000))
text(x = b, dat[1,] / 2, labels = dat[1,], cex = 0.8)
text(x = b, dat[1,] + dat[2,] - 100, labels = dat[2,], cex = 0.8)
text(x = b, colSums(dat) + 1000, labels = dat[3,], cex = 0.8)
legend('bottomright', c("both", 'none', "partial"), fill = c('red','blue', 'green'), 
       cex = 0.6, horiz = T)

#analysis of the kind of modification amongst ortholog pairs (no modification/complex_modif/1 domain modification)
dat = read.table(text = "BOVIN GORGO MACMU MONDO MOUSE PANTR PIGXX RATNO
1 1006  880   997   884  881   722   974   999
2 512   339   378   468  433   258   460   486
3 13087 13979 13043 8555 13676 14353 10214 12993", header = T)
b = barplot(as.matrix(dat), col = c('red', 'green', 'blue'),
            ylab = 'Pair count', main = 'Analysis of the kind of modification amongst ortholog pairs',
            cex.main = 0.7, cex.names = 0.8, ylim = c(0, max(dat) + 2000))
text(x = b, dat[1,] / 2, labels = dat[1,], cex = 0.8)
text(x = b, dat[1,] + dat[2,] + 300, labels = dat[2,], cex = 0.8)
text(x = b, colMeans(dat), labels = dat[3,], cex = 0.8)
legend('topright', c("1 domain modification", "complex modification", 'no modification/control'), fill = c('red','blue', 'green'), 
       cex = 0.6, horiz = F)

#analysis of domain repetition in 1 domain modification group (no repeat/domain repetition)
dat = read.table(text = "BOVIN GORGO MACMU MONDO MOUSE PANTR PIGXX RATNO
1 472 399 480 424 392 300 544 473
2 534 481 517 460 489 422 430 526", header = T)
b = barplot(as.matrix(dat), col = c('red', 'blue'),
            ylab = 'Pair count', main = 'Analysis of domain repetition in 1 domain modification group',
            cex.main = 0.7, cex.names = 0.8)
text(x = b, dat[1,] / 2, labels = dat[1,], cex = 0.8)
text(x = b, dat[1,] + dat[2,] / 2, labels = dat[2,], cex = 0.8)
legend('topright', c("domain loss not involved in repetition", "domain loss involved in repetition"), fill = c('red','blue'), 
       cex = 0.6, horiz = F)

#Ortholog with one domain modification in N- or C- termini called by domainDIFF
dat = read.table(text = "BOVIN GORGO MACMU MONDO MOUSE PANTR PIGXX RATNO
1 224786 232346 257597 183774 272425 206993 228663 242116", header = T)
b = barplot(as.matrix(dat), col = c('red'),
            ylab = 'Pair count', main = 'Ortholog with one domain modification in N- or C- termini called by domainDIFF',
            cex.main = 0.7, cex.names = 0.8)
text(x = b, dat[1,] / 2, labels = dat[], cex = 0.8)

#Intersect between domainDIFF and 1 domain not repeated modification group
dat = read.table(text = "BOVIN GORGO MACMU MONDO MOUSE PANTR PIGXX RATNO
1 352 319 374 319 276 242 426 354", header = T)
b = barplot(as.matrix(dat), col = c('red'),
            ylab = 'Pair count', main = 'Intersect between domainDIFF and 1 domain not repeated modification group',
            cex.main = 0.7, cex.names = 0.8)
text(x = b, dat[1,] / 2, labels = dat[], cex = 0.8)