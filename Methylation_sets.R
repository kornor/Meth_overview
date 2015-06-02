setwd("~/Bioinformatics Work/Meth & RNA/Meth_overview")

library(survival)



### Load in the clinical data
clin <- read.table("Clinical_final.txt", sep = "\t", header = TRUE, row.names = 1)

### Load in the exp data

exp <- read.table("Exp_final.txt", sep = "\t", header = TRUE, row.names = 1)

exp <- as.data.frame(t(exp))

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

### Seperate into low and high methylation groups by quartile
colnames(Fang_Ave)[1] <- "MethScore"
attach(Fang_Ave)
### for 4 quartiles

brks <- with(Fang_Ave, quantile(MethScore, probs = c(0, 0.25, 0.5, 0.75, 1)))
Values <- within(Fang_Ave, quartile <- cut(MethScore, breaks = brks, labels = 1:4, 
                                         include.lowest = TRUE))


library(plyr) ## for the "count" function, which is sweeeet!
count(Values$quartile)
## you can see that the distribution is very even / same kinda numbers in each quartile!

### for only high vs low

brks <- with(Fang_Ave, quantile(MethScore, probs = c(0, 0.5, 1)))
Values <- within(Fang_Ave, half <- cut(MethScore, breaks = brks, labels = 1:2, 
                                     include.lowest = TRUE))

count(Values$half)


#### and you're off and running!!

#### Set the MethScore quartiles as a factor

Quar <- as.factor(Values$quartile)

#### Look at expression of DNMT by quartile

boxplot(exp$DNMT1)



#######################################################################################
## load in the stir sets also

stir_beta <- read.table("Stir_beta.txt", sep = "\t", header = TRUE, row.names = 1)

stir_bin <- read.table("Stir_binary.txt", sep = "\t", header = TRUE, row.names = 1)

