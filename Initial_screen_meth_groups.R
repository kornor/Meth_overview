setwd("~/Bioinformatics Work/Meth & RNA/Meth_overview")

library(survival)


#############   Including normals first

### Load in the normals files for each subset

Fang_200a <- read.table("Fang_TSS200_all_exp_set.txt", sep ="\t", header = TRUE, row.names = 1)

########### create factors

Pam50 <- as.factor(Fang_200a$Pam50)
Quar <- as.factor(Fang_200a$Quartile)
Tert <- as.factor(Fang_200a$Tertile)


boxplot(Fang_200a$DNMT1 ~ Pam50)
boxplot(Fang_200a$DNMT1 ~ Quar)
boxplot(Fang_200a$DNMT1 ~ Tert)

boxplot(Fang_200a$MethScore ~ Pam50)

aov.out = aov(Fang_200a$DNMT3B ~ Pam50)
summary(aov.out)
TukeyHSD(aov.out)


### Survivial curves

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Quar),
                    data=Fang_200a) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Quar)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Quar,
               data=Fang_200a,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation quartile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Quar)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Quar), data = Fang_200a))


####### by tertile

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Tert),
                    data=Fang_200a) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Tert)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Tert,
               data=Fang_200a,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Terttile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Tert)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Tert), data = Fang_200a))


##########################################################################

Fang_200i <- read.table("Fang_TSS200_island_exp_set.txt", sep ="\t", header = TRUE, row.names = 1)

########### create factors

Pam50 <- as.factor(Fang_200i$Pam50)
Quar <- as.factor(Fang_200i$Quartile)
Tert <- as.factor(Fang_200i$Tertile)


boxplot(Fang_200i$DNMT1 ~ Pam50)
boxplot(Fang_200i$DNMT1 ~ Quar)
boxplot(Fang_200i$DNMT1 ~ Tert)

boxplot(Fang_200i$MethScore ~ Pam50)

aov.out = aov(Fang_200i$DNMT3B ~ Quar)
summary(aov.out)
TukeyHSD(aov.out)


### Survivial curves

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Quar),
                    data=Fang_200i) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Quar)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Quar,
               data=Fang_200i,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation quartile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Quar)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Quar), data = Fang_200i))


####### by tertile

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Tert),
                    data=Fang_200i) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Tert)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Tert,
               data=Fang_200i,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Terttile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Tert)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Tert), data = Fang_200i))

#################################################################################################

Fang_1500a <- read.table("Fang_TSS1500_all_exp_set.txt", sep ="\t", header = TRUE, row.names = 1)


########### create factors

Pam50 <- as.factor(Fang_1500a$Pam50)
Quar <- as.factor(Fang_1500a$Quartile)
Tert <- as.factor(Fang_1500a$Tertile)


boxplot(Fang_1500a$DNMT1 ~ Pam50)
boxplot(Fang_1500a$DNMT1 ~ Quar)
boxplot(Fang_1500a$DNMT1 ~ Tert)

boxplot(Fang_1500a$MethScore ~ Pam50)

aov.out = aov(Fang_1500a$DNMT3B ~ Quar)
summary(aov.out)
TukeyHSD(aov.out)


### Survivial curves

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Quar),
                    data=Fang_1500a) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Quar)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Quar,
               data=Fang_1500a,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation quartile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Quar)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Quar), data = Fang_1500a))


####### by tertile

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Tert),
                    data=Fang_1500a) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Tert)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Tert,
               data=Fang_1500a,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Terttile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Tert)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Tert), data = Fang_1500a))




#####################################################
Fang_1500i <- read.table("Fang_TSS1500_island_exp_set.txt", sep ="\t", header = TRUE, row.names = 1)



########### create factors

Pam50 <- as.factor(Fang_1500i$Pam50)
Quar <- as.factor(Fang_1500i$Quartile)
Tert <- as.factor(Fang_1500i$Tertile)


boxplot(Fang_1500i$DNMT1 ~ Pam50)
boxplot(Fang_1500i$DNMT1 ~ Quar)
boxplot(Fang_1500i$DNMT1 ~ Tert)

boxplot(Fang_1500i$MethScore ~ Pam50)

aov.out = aov(Fang_1500i$DNMT3B ~ Tert)
summary(aov.out)
TukeyHSD(aov.out)


### Survivial curves

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Quar),
                    data=Fang_1500i) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Quar)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Quar,
               data=Fang_1500i,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation quartile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Quar)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Quar), data = Fang_1500i))


####### by tertile

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Tert),
                    data=Fang_1500i) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Tert)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Tert,
               data=Fang_1500i,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Terttile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Tert)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Tert), data = Fang_1500i))


############################################
#### Stir groups


### Load in the normals files for each subset

Stir_200a <- read.table("Stir_TSS200_all_exp_set.txt", sep ="\t", header = TRUE, row.names = 1)

########### create factors

Pam50 <- as.factor(Stir_200a$Pam50)
Quar <- as.factor(Stir_200a$Quartile)
Tert <- as.factor(Stir_200a$Tertile)


boxplot(Stir_200a$DNMT1 ~ Pam50)
boxplot(Stir_200a$DNMT1 ~ Quar)
boxplot(Stir_200a$DNMT1 ~ Tert)

boxplot(Stir_200a$MethScore ~ Pam50)

aov.out = aov(Stir_200a$DNMT1 ~ Quar)
summary(aov.out)
TukeyHSD(aov.out)


### Survivial curves

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Quar),
                    data=Stir_200a) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Quar)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Quar,
               data=Stir_200a,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation quartile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Quar)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Quar), data = Stir_200a))


####### by tertile

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Tert),
                    data=Stir_200a) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Tert)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Tert,
               data=Stir_200a,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Terttile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Tert)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Tert), data = Stir_200a))


##########################################################################

Stir_200i <- read.table("Stir_TSS200_island_exp_set.txt", sep ="\t", header = TRUE, row.names = 1)

########### create factors

Pam50 <- as.factor(Stir_200i$Pam50)
Quar <- as.factor(Stir_200i$Quartile)
Tert <- as.factor(Stir_200i$Tertile)


boxplot(Stir_200i$DNMT1 ~ Pam50)
boxplot(Stir_200i$DNMT1 ~ Quar)
boxplot(Stir_200i$DNMT1 ~ Tert)

boxplot(Stir_200i$MethScore ~ Pam50)

aov.out = aov(Stir_200i$DNMT1 ~ Quar)
summary(aov.out)
TukeyHSD(aov.out)


### Survivial curves

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Quar),
                    data=Stir_200i) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Quar)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Quar,
               data=Stir_200i,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation quartile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Quar)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Quar), data = Stir_200i))


####### by tertile

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Tert),
                    data=Stir_200i) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Tert)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Tert,
               data=Stir_200i,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Terttile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Tert)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Tert), data = Stir_200i))

#################################################################################################

Stir_1500a <- read.table("Stir_TSS1500_all_exp_set.txt", sep ="\t", header = TRUE, row.names = 1)


########### create factors

Pam50 <- as.factor(Stir_1500a$Pam50)
Quar <- as.factor(Stir_1500a$Quartile)
Tert <- as.factor(Stir_1500a$Tertile)


boxplot(Stir_1500a$DNMT1 ~ Pam50)
boxplot(Stir_1500a$DNMT1 ~ Quar)
boxplot(Stir_1500a$DNMT1 ~ Tert)

boxplot(Stir_1500a$MethScore ~ Pam50)

aov.out = aov(Stir_1500a$DNMT1 ~ Quar)
summary(aov.out)
TukeyHSD(aov.out)


### Survivial curves

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Quar),
                    data=Stir_1500a) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Quar)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Quar,
               data=Stir_1500a,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation quartile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Quar)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Quar), data = Stir_1500a))


####### by tertile

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Tert),
                    data=Stir_1500a) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Tert)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Tert,
               data=Stir_1500a,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Terttile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Tert)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Tert), data = Stir_1500a))




#####################################################
Stir_1500i <- read.table("Stir_TSS1500_island_exp_set.txt", sep ="\t", header = TRUE, row.names = 1)



########### create factors

Pam50 <- as.factor(Stir_1500i$Pam50)
Quar <- as.factor(Stir_1500i$Quartile)
Tert <- as.factor(Stir_1500i$Tertile)


boxplot(Stir_1500i$DNMT1 ~ Pam50)
boxplot(Stir_1500i$DNMT1 ~ Quar)
boxplot(Stir_1500i$DNMT1 ~ Tert)

boxplot(Stir_1500i$MethScore ~ Pam50)

aov.out = aov(Stir_1500i$DNMT1 ~ Quar)
summary(aov.out)
TukeyHSD(aov.out)


### Survivial curves

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Quar),
                    data=Stir_1500i) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Quar)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Quar,
               data=Stir_1500i,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation quartile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Quar)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Quar), data = Stir_1500i))


####### by tertile

fit.diff = survdiff(Surv(Survival_time,Survival_event == 1) ~ factor(Tert),
                    data=Stir_1500i) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Tert)))-1),3) 
fit1 = survfit(Surv(Survival_time,Survival_event == 1)~Tert,
               data=Stir_1500i,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n methylation Terttile"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Tert)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(Survival_time,Survival_event == 1)~factor(Tert), data = Stir_1500i))


