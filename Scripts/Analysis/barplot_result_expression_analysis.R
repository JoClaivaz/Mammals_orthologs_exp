'''
Joaquim Claivaz
171003

plot rsults of expression number in mammals analysis
'''

##number of ortholog pair present in 1 domain modification or control group having expression data
#bash command on expression file parsed:
#cat Specie2_expression_parsed | cut -f3,7 | sort -u | grep control | wc -l
#cat Specie2_expression_parsed | cut -f3,7 | sort -u | grep modif | wc -l
dat = read.table(text = "BOVIN GORGO MACMU MONDO MOUSE PANTR PIGXX RATNO
1 13087 13979 13043 8555 13662 14353 10214 12993
2 352 319 374 319 276 242 426 354", header = T)
b = barplot(as.matrix(dat), col = c('blue', 'red'),
            ylab = 'Pair count', main = 'number of ortholog pair present in 1 domain modification or control group having expression data\nmissing data in parenthesis',
            cex.main = 0.7, cex.names = 0.8, ylim = c(0, 16000))
text(x = b, dat[1,] / 2, labels = dat[1,], cex = 0.8)
text(x = b, dat[1,] + dat[2,] + 1000, labels = dat[2,], cex = 0.8)
text(x = b, dat[1,] / 2 - 1000, labels = c(rep('', 4), '(14)', rep('', 3)), cex = 0.8)
legend('bottomright', c( "domain loss pair with expression data", "control pair with expression data"), fill = c('red','blue'), 
       cex = 0.6, horiz = F)
