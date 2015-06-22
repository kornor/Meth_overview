setwd("~/Bioinformatics Work/Meth & RNA/Meth_overview")
##############
#### This file cleans and preps the meth file inlcuding the normal samples
### Including the normals will result in changes to the designation of the tertile etc
### which might effect the results, hence keeping two seperate files


##### Read in exp & clin file
exp <- read.table("Exp_inc_normals.txt", sep = "\t", header = TRUE, row.names = 1)
exp <- as.data.frame(t(exp))

clin <- read.table("Clinical_final_normals.txt", sep = "\t", header = TRUE, row.names = 1)

### Use the meth file that already has duplicate removed
Meth1 <- read.table("Fang_TSS200_all.txt", sep = "\t", header = TRUE)

######## read in the normal file
Norm1<- read.table("Fang_TSS200_all_norm.txt", sep = "\t", header = TRUE, row.names = 1)

##### check the colnames are matching and remove excess
list <- intersect(colnames(Meth1), colnames(Norm1))
Meth3 <- Meth1[,list]


### bind together as a set
Meth4 <- rbind(Norm1, Meth3)

##### Match rownames to exp 
#### check samples match between files and trim

list <- intersect(row.names(Meth4), row.names(exp))
Meth4 <- Meth4[list,]
clin <- clin[list,]
exp <- exp[list,]

#### save out the new meth file

write.table(Meth4, "Fang_TSS200_all_norm_final.txt", sep = "\t")

######## set up the colmeans for meth set

# center with 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# apply it
CentreMeth<- center_colmeans(Meth4)

Meth_Ave <- as.data.frame(rowMeans(CentreMeth))

rownames(Meth_Ave) <- rownames(Meth4) 
colnames(Meth_Ave)[1] <- "MethScore"
#write.table(Meth_Ave, "Meth_beta_islands_ave.txt", sep = "\t")


### for 4 quartiles

brks <- with(Meth_Ave, quantile(MethScore, probs = c(0, 0.25, 0.5, 0.75, 1)))
Values <- within(Meth_Ave, quartile <- cut(MethScore, breaks = brks, labels = 1:4, 
                                           include.lowest = TRUE))


brks <- with(Meth_Ave, quantile(MethScore, probs = c(0, 0.33, 0.66, 1)))
Terts <- within(Meth_Ave, terts <- cut(MethScore, breaks = brks, labels = 1:3, 
                                       include.lowest = TRUE))

count(Values$quartile)
count(Terts$terts)


########### Set up dataframe for exports



look <- c("DNMT1", "DNMT3A", "DNMT3L", "DNMT3B", "UHRF1","MTHFR")
list <- match(look, names(exp))
exp_set <- exp[,list]

exp_set$Pam50 <- clin$PAM50Call_RNAseq
exp_set$Survival_time <- clin$OS_Time_nature2012
exp_set$Survival_event <- clin$OS_event_nature2012
exp_set$MethScore <- Meth_Ave$MethScore
exp_set$Quartile <- Values$quartile
exp_set$Tertile <- Terts$terts

write.table(exp_set, "Fang_TSS200_all_norm_exp_set.txt", sep = "\t")
#######################################################################################

### Use the meth file that already has duplicate removed
Meth1 <- read.table("Fang_TSS1500_all.txt", sep = "\t", header = TRUE)

######## read in the normal file
Norm1<- read.table("Fang_TSS1500_all_norm.txt", sep = "\t", header = TRUE, row.names = 1)

##### check the colnames are matching and remove excess
list <- intersect(colnames(Meth1), colnames(Norm1))
Meth3 <- Meth1[,list]

### bind together as a set
Meth4 <- rbind(Norm1, Meth3)

##### Match rownames to exp 
#### check samples match between files and trim

list <- intersect(row.names(Meth4), row.names(exp))
Meth4 <- Meth4[list,]
clin <- clin[list,]
exp <- exp[list,]

#### save out the new meth file

write.table(Meth4, "Fang_TSS1500_all_norm_final.txt", sep = "\t")

######## set up the colmeans for meth set

# center with 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# apply it
CentreMeth<- center_colmeans(Meth4)

Meth_Ave <- as.data.frame(rowMeans(CentreMeth))

rownames(Meth_Ave) <- rownames(Meth4) 
colnames(Meth_Ave)[1] <- "MethScore"
#write.table(Meth_Ave, "Meth_beta_islands_ave.txt", sep = "\t")


### for 4 quartiles

brks <- with(Meth_Ave, quantile(MethScore, probs = c(0, 0.25, 0.5, 0.75, 1)))
Values <- within(Meth_Ave, quartile <- cut(MethScore, breaks = brks, labels = 1:4, 
                                           include.lowest = TRUE))


brks <- with(Meth_Ave, quantile(MethScore, probs = c(0, 0.33, 0.66, 1)))
Terts <- within(Meth_Ave, terts <- cut(MethScore, breaks = brks, labels = 1:3, 
                                       include.lowest = TRUE))

count(Values$quartile)
count(Terts$terts)


########### Set up dataframe for exports



look <- c("DNMT1", "DNMT3A", "DNMT3L", "DNMT3B", "UHRF1","MTHFR")
list <- match(look, names(exp))
exp_set <- exp[,list]

exp_set$Pam50 <- clin$PAM50Call_RNAseq
exp_set$Survival_time <- clin$OS_Time_nature2012
exp_set$Survival_event <- clin$OS_event_nature2012
exp_set$MethScore <- Meth_Ave$MethScore
exp_set$Quartile <- Values$quartile
exp_set$Tertile <- Terts$terts

write.table(exp_set, "Fang_TSS1500_all_norm_exp_set.txt", sep = "\t")
#######################################################################################
