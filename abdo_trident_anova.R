#----------------------------------------------------------------------
#author marlies dolezal
#november 2016
# changed: Lukas Endler
# Tues. Nov. 14th 2017
#----------------------------------------------------------------------


setwd("~/Data/Trident")

getwd()
#file<-read.table("tan-cuticles-for-Openstat.red2.csv", header = T, skip=0, sep = ";",na.strings=".")
library(lsmeans)
library("multcomp")
#install.packages("multcompView")
library("multcompView")

library(car)
library(MASS)
#install.packages("estimability")
library("estimability")

#function to calculate eta squared and partial eta squared
ETA2 <- function(lm.obj){
  #cat("\n       ETA2:\n\n")
  tmp <- Anova(lm.obj,type=3)
  print("Total Sum Sq")
  print(sum(tmp$'Sum Sq'[2:9]))
  ETA2<-round(tmp$'Sum Sq'[2:9]/sum(tmp$'Sum Sq'[2:9]),4)
  partialETA2<-round(tmp$'Sum Sq'[2:9]/(tmp$'Sum Sq'[2:9] + tmp$'Sum Sq'[9]),4)
  SumSqlabel<-rownames(tmp)[2:9]
  print(rbind(SumSqlabel,ETA2,partialETA2))
  return(data.frame("SumSqlabel" = SumSqlabel,"ETA2"= ETA2,"pETA2"=partialETA2))
}

#----------------------

#head(dl[[dn]])
abdo<-read.table("tan-cuticles-for-Openstat.csv", header = T, skip=0, sep = ";")
trident<-read.table("Pigmentation-trident.csv", header = T, skip=0, sep = ";")
names(trident)[5]="trident"
dl=list(abdo,trident)
names(dl)=c("abdo","trident")

for (di in 1:length(dl)){
  for (field in c("SNP1","SNP2","SNP3")){
    dl[[di]][,field]<-as.factor(ifelse(dl[[di]][,field]==1,"D","L"))
    #dl[[di]][,field]<-relevel(dl[[di]][,field],"L")
  }
}


abdo=dl[["abdo"]]
trident=dl[["trident"]]

#detach(abdo)
#loop through segment A5-A7 as phenotypes
contrasts(dl[["abdo"]]$SNP1) = contr.sum(2)
contrasts(dl[["abdo"]]$SNP2) = contr.sum(2)
contrasts(dl[["abdo"]]$SNP3) = contr.sum(2)
contrasts(dl[["trident"]]$SNP1) = contr.sum(2)
contrasts(dl[["trident"]]$SNP2) = contr.sum(2)
contrasts(dl[["trident"]]$SNP3) = contr.sum(2)

pig.lm1=list()
pig.lm2=list()

etas.df <- NULL

for (dn in names(dl)) {
  for (i in c(names(abdo)[6:8],"trident")) {
    if ( ! i %in% names(dl[[dn]]) ) {next}  
    print("                                           ")
    print("###########################################")
    print(dl[[dn]][,i])
    print("###########################################")
    
    print("----------------------------------------------------")
    print("MAIN EFFECTS / ADDITIVE MODEL - sum zero restriction")
    print("----------------------------------------------------")
    pig.lm1[[i]]<-lm(dl[[dn]][,i]~SNP1+SNP2+SNP3, data=dl[[dn]], contrasts="contr.sum")
    print(summary(pig.lm1[[i]]))
    print(Anova(pig.lm1[[i]],type=3))
    #print(model.matrix(lm1))
    #print(contrasts(abdo$SNP_1))
    #print(contrasts(abdo$SNP_2))
    #print(contrasts(abdo$SNP_3))
    
    print("------------------------------------------------------------------")
    print("FULL FACTORIAL MODEL / ADDITIVE & EPISTATIC - sum zero restriction")
    print("------------------------------------------------------------------")
    
    pig.lm2[[i]]<-lm(dl[[dn]][,i]~SNP1*SNP2*SNP3, data=dl[[dn]], contrasts="contr.sum")
    #print(model.matrix(lm2))
    print(summary(pig.lm2[[i]]))
    print("ANOVA type 3 SNP1*SNP2*SNP3  contr.sum")
    print(Anova(pig.lm2[[i]],type=3))
    #print(model.matrix(lm1))
    
    print("----------------------------------------------------")
    print("comparing main effects and full factorial model     ")
    print("----------------------------------------------------")
    print(anova(pig.lm1[[i]],pig.lm2[[i]]))
    
    print("----------------------------------------------------")
    print("ETA^2 and partial ETA^2 full factorial model     ")
    print("----------------------------------------------------")
    etas <- ETA2(pig.lm2[[i]])
    etas$P <- round(Anova(pig.lm2[[i]],type=3)$"Pr(>F)"[2:9],3)
    colnames(etas) <- sapply(colnames(etas),paste,i,sep="_")
    if (is.null(etas.df)){
      etas.df <- etas[,2:4]
      rownames(etas.df) <- etas[,1]
    } else {
      etas.df <- cbind(etas.df,etas[,2:4])
    }
  }
}


qqPlot(pig.lm2[[i]])
lsmeans::cld(lsmeans(pig.lm2[[i]],c("SNP1","SNP2","SNP3")))
lsmip(pig.lm2[[i]],SNP1~SNP2|SNP3)
plot(lsmeans(pig.lm2[[i]],~SNP1))
ref.grid(pig.lm2[[i]])
summary(ref.grid(pig.lm2[[i]]))
with(trident, tapply(trident, SNP1, mean))
with(trident, tapply(trident, SNP2, mean))
with(trident, tapply(trident, SNP3, mean))
lsmeans(pig.lm2[[i]],~ SNP1)
contrast(lsmeans(pig.lm2[[i]],~ SNP1),method="pairwise")
confint(contrast(lsmeans(pig.lm2[[i]],~ SNP1),method="pairwise"))
Anova(pig.lm2[[i]],type="III")
leveneTest(pig.lm2[[i]])

labs=c()
for (i in rev(levels(dl[["trident"]]$SNP1))) { 
  for (j in rev(levels(dl[["trident"]]$SNP2))) { 
    for (l in rev(levels(dl[["trident"]]$SNP3))) { 
      labs=c(labs,paste0(i,j,l)) }}}
labs=apply(expand.grid( levels(dl[["trident"]]$SNP1), levels(dl[["trident"]]$SNP2), levels(dl[["trident"]]$SNP3)  ),1, paste0, collapse = "")
# create boxplot
pdf("boxplot_trident_pig.pdf",width=7, height=5)
boxplot( trident ~ - SNP1 - SNP2 - SNP3 ,dl[["trident"]],names=labs)
#boxplot( trident ~ + SNP1 + SNP2 + SNP3 ,dl[["trident"]])
dev.off()



