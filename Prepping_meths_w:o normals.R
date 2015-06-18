setwd("~/Bioinformatics Work/Meth & RNA/Meth_overview")

#########
## This file takes the meth data and cleans and preps it for later analysis in prism
####  After completion of running this file - find & replace "all" with "island"
#### Then find and replace "Fang" with "Stir" and do all over again.
#### Then proceed to do the normals combine file set

##### Read in exp & clin file
exp <- read.table("Exp_inc_normals.txt", sep = "\t", header = TRUE, row.names = 1)
exp <- as.data.frame(t(exp))

clin <- read.table("Clinical_final_normals.txt", sep = "\t", header = TRUE, row.names = 1)

################ read in the meth file, wihtout row names

Meth1 <- read.table("Fang_TSS200_all.txt", sep = "\t", header = TRUE)

### find dupliucated rownames and remove
list <- which(duplicated(Meth1[,1]))

Meth2 <- Meth2[-list,]

####### write out the table / open, fix in excel and reopen with rownames
write.table(Meth2, "Fang_TSS200_all.txt", sep = "\t")
Meth3 <- read.table("Fang_TSS200_all.txt", sep = "\t", header = TRUE, row.names = 1)

##### Match rownames to exp 
#### check samples match between files and trim

list <- intersect(row.names(Meth3), row.names(exp))
Meth3 <- Meth3[list,]

#### save out the new meth file

write.table(Meth3, "Fang_TSS200_all_final.txt", sep = "\t")

######## set up the colmeans for meth set

# center with 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# apply it
CentreMeth<- center_colmeans(Meth3)

Meth_Ave <- as.data.frame(rowMeans(CentreMeth))

rownames(Meth_Ave) <- rownames(Meth3) 
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

write.table(exp_set, "Fang_TSS200_all_exp_set.txt", sep = "\t")
#######################################################################################
###########  TSS1500 now

################ read in the meth file, wihtout row names

Meth1 <- read.table("Fang_TSS1500_all.txt", sep = "\t", header = TRUE)

### find dupliucated rownames and remove
list <- which(duplicated(Meth1[,1]))

Meth2 <- Meth2[-list,]

####### write out the table / open, fix in excel and reopen with rownames
write.table(Meth2, "Fang_TSS1500_all.txt", sep = "\t")
Meth3 <- read.table("Fang_TSS1500_all.txt", sep = "\t", header = TRUE, row.names = 1)

##### Match rownames to exp 
#### check samples match between files and trim

list <- intersect(row.names(Meth3), row.names(exp))
Meth3 <- Meth3[list,]

#### save out the new meth file

write.table(Meth3, "Fang_TSS1500_all_final.txt", sep = "\t")

######## set up the colmeans for meth set

# center with 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# apply it
CentreMeth<- center_colmeans(Meth3)

Meth_Ave <- as.data.frame(rowMeans(CentreMeth))

rownames(Meth_Ave) <- rownames(Meth3) 
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
  

write.table(exp_set, "Fang_TSS1500_all_exp_set.txt", sep = "\t")