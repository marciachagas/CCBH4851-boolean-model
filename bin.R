#Package installation of "Binarize"
install.packages("Binarize")

#Package loading
library("Binarize")

setwd("D:/Mestrado/Binarização/Genes")

dati_csv <-read.csv("D:/Mestrado/Binarização/Genes/Genes-separados/ansR.csv",header=T)

dati_matrix <- as.matrix(t(dati_csv))

binMatrix <- binarizeMatrix(dati_matrix,
                            method="BASCB",
                            adjustment="fdr",
                            tau=0.5)


print(binMatrix)


write.csv(binMatrix ,file="D:/Mestrado/Binarização/Genes/Genes-separados/Genes-bin/ansR-bin.csv")

