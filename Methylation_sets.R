setwd("~/Bioinformatics Work/Meth & RNA/Meth_overview")

library(survival)



### Load in the clinical data
clin <- read.table("Clinical_final_normals.txt", sep = "\t", header = TRUE, row.names = 1)

### Load in the exp data

expt <- read.table("Exp_final.txt", sep = "\t", header = TRUE, row.names = 1)

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
colnames(Fang_Ave)[1] <- "MethScore"
write.table(Fang_Ave, "Fang_beta_average.txt", sep = "\t")

hist(Fang_Ave$`rowMeans(CentreFang)`)

#########################################

##Do averages of the binary

Fang_Collect <- as.data.frame(rowMeans(fang_bin))
colnames(Fang_Collect)[1] <- "MethBin"
hist(Fang_Collect$MethBin)


##########################################

## Only the top 3 Fang CIMP genes
B_cimp <- data.frame(fang_bin$ALX4)
B_cimp$ARHGEF7 <- fang_bin$ARHGEF7
B_cimp$RASGRF2 <- fang_bin$RASGRF2

colnames(B_cimp)[1] <- "ALX4"

### Get a sum for the 3 genes methylation status
B_cimp$Total <- rowSums(B_cimp) 
## Set the cut off (2/3 genes methylated = B_CIMP+)
B_cimp$ME_Rank[B_cimp$Total>=2] <-"Pos"
B_cimp$ME_Rank[B_cimp$Total<=1]<- "Neg"

count(B_cimp$ME_Rank)
##Most are neg??

##set as a factor

CIMP <- as.factor(B_cimp$ME_Rank)

### Seperate into low and high methylation groups by quartile - beta group

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
Half <- as.factor(Values$half)

#### Look at expression of DNMT by quartile
### Do a box plot, anova, etc

boxplot(exp$DNMT1 ~ Quar)

### AOV testing shows that DNMT is higher in basal  than other subtype
aov.out = aov(exp$DNMT1 ~ Pam50)
summary(aov.out)
TukeyHSD(aov.out)


### And only marginally higher in the lowest quartile of FANG methylation (p = ~ 0.1)
aov.out = aov(exp$DNMT1 ~ Quar)
summary(aov.out)
TukeyHSD(aov.out)
##CIMP is not related to DNMT expression
aov.out = aov(exp$DNMT1 ~ CIMP)
summary(aov.out)
TukeyHSD(aov.out)


### What about methylation quartile by Pam50?

boxplot(Fang_Ave$MethScore ~ Pam50)
aov.out = aov(Fang_Ave$MethScore ~ Pam50)
summary(aov.out)
TukeyHSD(aov.out)
### This shows that methylation is reduced in basals com to others 
## What about in Stir?
boxplot(Stir_Ave$MethScore ~ Pam50)
aov.out = aov(Stir_Ave$MethScore ~ Pam50)
summary(aov.out)
TukeyHSD(aov.out)
### Pattern is still there but not as pronounced (as expected). 


#### DNMT1 levels are more strongly tied to the basal phenotype than to the methylation pattern???



boxplot(exp$DNMT1 ~ Half)
boxplot(exp$DNMT3A ~ CIMP)

boxplot(exp$DNMT3A ~ Quar)

boxplot(exp$DNMT3B ~ Quar)

boxplot(exp$AHCY ~ Quar)

boxplot(exp$MAT1A ~ Quar)

### Survivial curves

fit.diff = survdiff(Surv(OS_Time_nature2012,OS_event_nature2012 == 1) ~ factor(Pam50),
                    data=clin) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Pam50)))-1),3) 
fit1 = survfit(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~Pam50,
               data=clin,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Pam50"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Pam50)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~factor(CIMP), data = clin))




#######################################################################################
## load in the stir sets also

stir_beta <- read.table("Stir_beta.txt", sep = "\t", header = TRUE, row.names = 1)

stir_bin <- read.table("Stir_binary.txt", sep = "\t", header = TRUE, row.names = 1)

### Do mean centred averages of the beta ############################

# center with 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# apply it
CentreStir <- center_colmeans(stir_beta)

Stir_Ave <- as.data.frame(rowMeans(CentreStir))

rownames(Stir_Ave) <- rownames(stir_beta) 
colnames(Stir_Ave)[1] <- "MethScore"

write.table(Stir_Ave, "Stir_beta_average.txt", sep = "\t")

hist(Stir_Ave$MethScore)

#########################################

##Do averages of the binary

Stir_Collect <- as.data.frame(rowMeans(stir_bin))
colnames(Stir_Collect)[1] <- "MethBin"

hist(Stir_Collect$MethBin)


#######################
attach(Stir_Ave)
### for 4 quartiles

brks <- with(Stir_Ave, quantile(MethScore, probs = c(0, 0.25, 0.5, 0.75, 1)))
Values <- within(Stir_Ave, quartile <- cut(MethScore, breaks = brks, labels = 1:4, 
                                           include.lowest = TRUE))


library(plyr) ## for the "count" function, which is sweeeet!
count(Values$quartile)
## you can see that the distribution is very even / same kinda numbers in each quartile!

### for only high vs low

brks <- with(Stir_Ave, quantile(MethScore, probs = c(0, 0.5, 1)))
Values <- within(Stir_Ave, half <- cut(MethScore, breaks = brks, labels = 1:2, 
                                       include.lowest = TRUE))

count(Values$half)


#### and you're off and running!!

#### Set the MethScore quartiles as a factor

Quar <- as.factor(Values$quartile)
Half <- as.factor(Values$half)

boxplot(exp$DNMT1 ~ Half)

### Survivial curves

fit.diff = survdiff(Surv(OS_Time_nature2012,OS_event_nature2012 == 1) ~ factor(Half),
                    data=clin) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Half)))-1),3) 
fit1 = survfit(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~Half,
               data=clin,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Halftile "), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Half)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~factor(Quar), data = clin))




############## Methylation by subtype

Pam50 <- factor(clin$PAM50Call_RNAseq) 
ER <- factor(clin$ER_Status_nature2012)
PR <- factor(clin$PR_Status_nature2012)
Her <- factor(clin$HER2_Final_Status_nature2012)
Mets <- factor(clin$Metastasis_Coded_nature2012)
Node <- factor(clin$Node_Coded_nature2012)


stripchart(Fang_Ave$MethScore ~ Pam50)
stripchart(Stir_Ave$MethScore ~ Pam50)

stripchart(Fang_Collect$MethBin ~ Pam50)
stripchart(Stir_Collect$MethBin ~ Pam50)

stripchart(Fang_Ave$MethScore ~ ER)
stripchart(Stir_Ave$MethScore ~ Her)


boxplot(exp$DNMT1 ~ Pam50)
boxplot(Fang_Ave$MethScore ~ Pam50)
boxplot(Stir_Ave$MethScore ~ Pam50)



############# What about the predictive value just within the basal subtype?
clinbasal <- subset(clin, Pam50 == "Basal")
stirbasal <- subset(Stir_Ave, Pam50 == "Basal")
fangbasal <- subset(Fang_Ave, Pam50 == "Basal")

brks <- with(stirbasal, quantile(MethScore, probs = c(0, 0.5, 1)))
stirbasal <- within(stirbasal, half <- cut(MethScore, breaks = brks, labels = 1:2, 
                                       include.lowest = TRUE))

brks <- with(fangbasal, quantile(MethScore, probs = c(0, 0.5, 1)))
fangbasal <- within(fangbasal, half <- cut(MethScore, breaks = brks, labels = 1:2, 
                                                include.lowest = TRUE))


brks <- with(stirbasal, quantile(MethScore, probs = c(0, 0.25, 0.5, 0.75, 1)))
stirbasal <- within(stirbasal, quartile <- cut(MethScore, breaks = brks, labels = 1:4, 
                                           include.lowest = TRUE))

brks <- with(fangbasal, quantile(MethScore, probs = c(0, 0.25, 0.5, 0.75, 1)))
fangbasal <- within(fangbasal, quartile <- cut(MethScore, breaks = brks, labels = 1:4, 
                                              include.lowest = TRUE))

count(fangbasal$half)

count(stirbasal$half)

fangHalf <- as.factor(fangbasal$half)
stirHalf <- as.factor(stirbasal$half)

stirQuar <- as.factor(stirbasal$quartile)
fangQuar <- as.factor(fangbasal$quartile)


### Survivial curves

fit.diff = survdiff(Surv(OS_Time_nature2012,OS_event_nature2012 == 1) ~ factor(fangHalf),
                    data=clinbasal) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(fangHalf)))-1),3) 
fit1 = survfit(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~fangHalf,
               data=clinbasal,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation group "), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(fangHalf)),
       fill = kmcolours, cex = 1)

