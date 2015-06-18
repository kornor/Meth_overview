

################ read in the meth file, wihtout row names

fang_island <- read.table("Fang_Island.txt", sep = "\t", header = TRUE)

### find dupliucated rownames and remove
list <- which(duplicated(fang_island[,1]))

fang_island <- fang_island[-list,]

####### write out the table / open, fix in excel and reopen with rownames
write.table(fang_island, "Fang_Island.txt", sep = "\t")
fang_island <- read.table("Fang_Island.txt", sep = "\t", header = TRUE, row.names = 1)

######## read in the normal file
norm_islands <- read.table("Fang_normals_island.txt", sep = "\t", header = TRUE, row.names = 1)

##### check the colnames are matching and remove excess
list <- intersect(colnames(fang_island), colnames(norm_islands))
fang_island <- fang_island[,list]

### bind together as a set
meth <- rbind(norm_islands, fang_island)


#### read in clin if not already done
clin<- read.table("Clinical_normals_inc.txt", sep = "\t", header = TRUE, row.names = 1)

#### check samples match between files and trim

list <- intersect(row.names(meth), row.names(exp))
meth <- meth[list,]

clin <- clin[list,]

exp <- exp[list,]
exp <- as.data.frame(t(exp))

#### save out the new meth file

write.table(meth, "Fang_islands_inc_normal.txt", sep = "\t")

write.table(clin, "Clinical_normals_inc.txt", sep = "\t")

write.table(exp, "Exp_inc_normals.txt", sep = "\t")

exp <- read.table("Exp_inc_normals.txt", sep = "\t", header = TRUE, row.names =1)
exp <- t(exp)



######## set up the colmeans for normal-meth set

# center with 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# apply it
CentreMeth<- center_colmeans(meth)

Meth_Ave <- as.data.frame(rowMeans(CentreMeth))

rownames(Meth_Ave) <- rownames(meth) 
colnames(Meth_Ave)[1] <- "MethScore"
write.table(Meth_Ave, "Meth_beta_islands_ave.txt", sep = "\t")



hist(Meth_Ave$MethScore)


### for 4 quartiles

brks <- with(Meth_Ave, quantile(MethScore, probs = c(0, 0.25, 0.5, 0.75, 1)))
Values <- within(Meth_Ave, quartile <- cut(MethScore, breaks = brks, labels = 1:4, 
                                           include.lowest = TRUE))


brks <- with(Meth_Ave, quantile(MethScore, probs = c(0, 0.33, 0.66, 1)))
Terts <- within(Meth_Ave, terts <- cut(MethScore, breaks = brks, labels = 1:3, 
                                           include.lowest = TRUE))

count(Values$quartile)
count(Terts$terts)


library(plyr) ## for the "count" function, which is sweeeet!
count(Values$quartile)

#### Set the MethScore quartiles as a factor

Quar <- as.factor(Values$quartile)
Tert <- as.factor(Terts$terts)


Pam50 <- factor(clin$PAM50Call_RNAseq) 
ER <- factor(clin$ER_Status_nature2012)
PR <- factor(clin$PR_Status_nature2012)
Her <- factor(clin$HER2_Final_Status_nature2012)
Mets <- factor(clin$Metastasis_Coded_nature2012)
Node <- factor(clin$Node_Coded_nature2012)

#### DNMT levels by group

boxplot(exp$DNMT1 ~ Pam50)
boxplot(exp$DNMT3A ~ Pam50)
boxplot(exp$DNMT3L ~ Pam50)
boxplot(exp$DNMT3B ~ Pam50)

boxplot(exp$CEBPA ~ Pam50)
### DNMT levels by methylation - don't hold when the normals are included
#### probably becuase the normals are hypomethylation

boxplot(exp$DNMT1 ~ Tert)
boxplot(exp$DNMT3A ~ Quar)
boxplot(exp$CEBPA ~ Quar)

boxplot(exp$MTHFR ~ Pam50)

### Methylation by group

boxplot(Meth_Ave$MethScore ~ Pam50)



look <- c("DNMT1", "DNMT3A", "DNMT3L", "DNMT3B", "UHRF1","MTHFR")
list <- match(look, names(exp))
exp_set <- exp[,list]

exp_set$Pam50 <- clin$PAM50Call_RNAseq

############# Need to finalise the methylation values and add new columns (MethScore & Values)
############# Then import into prism. 
