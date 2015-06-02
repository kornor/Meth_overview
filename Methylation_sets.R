setwd("~/Bioinformatics Work/Meth & RNA/Meth_overview")

library(survival)



### Load in the clinical data
clin <- read.table(clin)



### Load in the Fang set as a beta

fang_beta <- read.table("Fang_beta.txt", sep = "\t", header = TRUE, row.names = 1)

## and as a binary

fang_bin <- read.table("fang_binary.txt", sep = "\t", header = TRUE, row.names = 1)

### Do mean centred averages of the beta ############################

# center with 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# apply it
CentreFang <- center_colmeans(fang_beta)

Fang_Ave <- as.data.frame(rowMeans(CentreFang))

rownames(Fang_Ave) <- rownames(fang_beta) 
write.table(Fang_Ave, "Fang_beta_average.txt", sep = "\t")

hist(Fang_Ave$`rowMeans(CentreFang)`)

#########################################

##Do averages of the binary

Fang_Collect <- as.data.frame(rowMeans(fang_bin))

hist(Fang_Collect$`rowMeans(fang_bin)`)


##########################################


