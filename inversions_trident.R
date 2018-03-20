library(doBy)

boxplot_inv <- function(a,n,tit="",col1="lightblue",col2="orange",name1="CS",name2="B")
{
	pops=(length(a[[n]])-5)/2
	names=1:(length(a[[n]])-5)
	for(i in 1:pops){names[i]=paste(name1,i,sep=""); names[i+pops]=paste(name2,i,sep="");}
	boxplot(a[[n]][6:length(a[[n]])],data=a[[n]],notch=F,col=c(rep(col1,pops),rep(col2,pops)),ylab="frequency of allele fixed in inversion",names=names, main=paste("Inversion",a[[n]]$CHR[1],a[[n]]$INV[1],tit))
	}

boxplot_inv_combine <- function(a,tit="",col1="lightblue",col2="orange",name1="CS",name2="B",pops1=1:3,pops2=4:6)
{	
	b1=list()
	b2=list()
	group_names = c()
	pops1=pops1+5
	pops2=pops2+5
	for(i in names(a)){
		b1[[i]]=unlist(a[[i]][,pops1],use.names=F)
		b2[[i]]=unlist(a[[i]][,pops2],use.names=F)
		group_names = append(group_names,paste("IN","(",a[[i]]$CHR[1],")",toupper
(a[[i]]$INV[1]),sep ="")) 
		}
	group_names = c( "IN(2L)t"  ,"IN(2R)NS" ,"IN(3L)P" , "IN(3R)C" , "IN(3R)Mo" ,"IN(3R)P")
	num_all = length(names(a))
	boxplot(b1,at = 0:(num_all-1)*3 + 1, xlim = c(0,num_all*3), ylim = range(b1, b2), xaxt = "n", ylab="frequency of allele fixed in inversion",col=col1, main=tit,notch=T)
    boxplot(b2, at = 0:(num_all-1)*3 + 2, xaxt = "n",col=col2, add = TRUE,notch=T)
   legend("topright",c(name1,name2),fill=c(col1,col2))
  axis(1, at = 0:(num_all-1)*3 + 1.5, labels = group_names, tick = TRUE)
	}


# cs inversion plotting
setwd("/Volumes/vetgrid10//Data/Trident/BGI_111_112/Joined")
#setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH")
inv_data_2011 <- read.table("CMH/males_join_pTLV1_pTLV2_pTDV1_pTDV2_pTLI1_pTLI2_pTLI3_pTDI1_pTDI2_pTDI3__q20_filt_mc2_chroms_only_mct8_mcv15.inversions",header=TRUE)
summaryBy(pop1 + pop2 + pop3 + pop4 +pop5 + pop6 + pop7 + pop8 + pop9 + pop10 ~ CHR + INV, data = inv_data_2011, FUN =c(median), keep.names = T)
#aggregate(cbind(pop1,pop2,pop3,pop4)~CHR+INV,data=inv_data,FUN=median)
head(inv_data_2011)
boxplot_inv(a_11,1,tit=", all_pops")
#setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_Poly/CMH")
pops=c("LV1","LV2","DV1","DV2","LI1","LI2","LI3","DI1","DI2","DI3")
df=NULL
for(i in seq(6,24,by=2)){
  subtb <- inv_data_2011[,c(1,2,3,4,5,i,i+1)]
  colnames(subtb) <- c(colnames(inv_data_2011)[1:5],"freq","cov")
  subtb$inv_count <- round(subtb$freq*subtb$cov,digits=0)
  subtb$not_inv_count<- subtb$cov-subtb$inv_count
  facs=strsplit(pops[(i-4)/2],"")[[1]]
  print(facs)
  subtb$pig=facs[1]
  subtb$inv=apply(inv_data_2011[,c(1,2)],1,paste0,collapse="|")
  subtb$locus=apply(inv_data_2011[,c(1,2,3)],1,paste0,collapse="_")
  subtb$pop=facs[2]
  subtb$repl=paste0(facs[c(2,3)],collapse="")
  subtb$rep=facs[3]
  subtb$exp=paste0(facs[c(1,2,3)],collapse="")
  if (is.null(df)) {df <- subtb} else { df=rbind(df,subtb) }
}
df$pig=factor(df$pig,levels=c("L","D"))
df$locus=factor(df$locus)
df$pop=factor(df$pop,levels=c("V","I"))
df$repl=factor(df$repl,levels=c("V1","V2","I1","I2","I3"))
df$rep=factor(df$rep,levels=c("1","2","3"))
df$exp=factor(df$exp)
df$inv=factor(df$inv)
summary(df)  
a_trident <- splitBy(formula=~ CHR + INV, data = df )

head(df)
med_trident = aggregate(freq~inv+pig+pop+rep,data = df[,c("freq","pig","inv","pop","rep")],median)
str(med_trident)
library("lme4")
library("lmerTest")
library("lattice")
med_trident$repl=paste0(med_trident$pop,med_trident$rep)

p.vals <- c()
p.vals.aov <- c()
lmer_med_trident=list()
for (inversion in levels(med_trident$inv) ) {
  boxplot( freq *100 ~ pig + repl   , data=a_trident[[inversion]],col=(c("yellow","brown")), main=inversion)
  
  #a <- lmer( freq ~ 1 + pig + pop + (1|repl), data = subset(med_trident,inv == inversion))
  a <- lm( freq ~ 1 + pop + pig, data = subset(med_trident,inv == inversion))
  print(inversion)
  print(summary(a))
  p.vals <- c(p.vals,summary(a)$coefficients["pigD","Pr(>|t|)"])
  names(p.vals)[length(p.vals)] <- inversion
  p.vals.aov <- c(p.vals,summary(aov(a))$coefficients["pigD","Pr(>F)"])
  
  #qqnorm(ranef(a)[[1]][,1],main=paste(c("QQ Plot for random intercept locus",inversion),collapse=" "))
  #qqline(ranef(a)[[1]][,1],col="red")
  #qqnorm(ranef(a)[[2]][,1],main=paste(c("QQ Plot for random intercept replicates",inversion),collapse=" "))
  #qqline(ranef(a)[[2]][,1],col="red")
  qqnorm(resid(a),main=paste(c("QQ Plot for residuals",inversion),collapse=" "))
  qqline(resid(a),col="red")
  plot(a,which=1)
  plot(a,which=2)
} 

p.vals
p.vals.aov
p.adjust(p.vals,method="bonferroni")

library("nlme")
p.vals <- c()
for (inversion in levels(med_trident$inv) ) {
  #boxplot( freq *100 ~ pig + repl   , data=a_trident[[inversion]],col=(c("yellow","brown")), main=inversion)
  #boxplot( freq *100 ~ pig  + pop  , data=subset(med_trident,inv == inversion),col=(c("yellow","brown")), main=inversion)
  # a <- lmer( freq ~ 1 + (1|pop) + (1|repl) + pig, data = subset(med_trident,inv == inversion))
  a <- lme( freq ~ 1 + pig, data = subset(med_trident,inv == inversion), random = ~ 1|pop)
  print(inversion)
  print(anova(a))
  p.vals <- c(p.vals,anova(a)[["p-value"]][2])
  names(p.vals)[length(p.vals)] <- inversion
  #print(summary(a))
} 
p.adjust(p.vals,method="BY")

arcsine <- function(frac) { asin(sqrt(frac))}
inv.arcsine <- function(p) { sin(p)^2}
logit <- function(frac) { log((frac+0.000000000001)/(1-frac)) }
inv.logit <- function(x){exp(x)/(1+exp(x))}

inversion="3L|p"
inversion="3R|p"
inversion="2L|t"
inversion="3R|c"
inversion="3R|mo"
inversion="2R|ns"

p.vals <- c()
p.vals_v <- c()
p.vals_i <- c()
a <- lmer( arcsine(freq) ~ 1 + pop + pig + (1|repl) + (1|locus), data = a_trident[[inversion]])
a <- lmer( arcsine(freq) ~ 1 + pop*pig + (1|repl) + (1|locus), data = a_trident[[inversion]])

a <- lmer( arcsine(freq) ~ 1 + pig + (1|repl) + (1|locus), data = subset(a_trident[[inversion]], pop == "V"))
library(reshape2)
dcast(a_trident[[inversion]],locus)

summary(a)
for (inversion in names(a_trident) ) {
  a <- lmer( arcsine(freq) ~ 1 + pop + pig + (1|repl) + (1|locus), data = a_trident[[inversion]])
  #a_nr <- lmer( arcsine(freq) ~ 0 + (1|pop/rep) + pig, data = a_trident[[inversion]])
  #a1 <- lmer( arcsine(freq) ~ 0 + (1|pop/rep) + (0+pig|pop/rep) + pig, data = a_trident[[inversion]])
  #a <- lmer( arcsine(freq) ~ 1 + (1 + pig |pop/rep)  + pig, data = a_trident[[inversion]])
  #p.vals=c(p.vals,summary(a)[["coefficients"]][2,5])
  #names(p.vals)[length(p.vals)]=inversion
  
  print(inversion)
 print(warnings())
 #print(anova(a_nr,a,refit=FALSE))
 randomeff <-ranef(a)
 # independence of random intercept & random slope
 # plot(randomeff[[1]][,1],randomeff[[2]][,1],
 #      xlab ="random intercept" , ylab ="random slope",
 #      main=paste(c("correlated intercept & slope model",inversion),collapse=" "))
 # # independence of random intercepts & residuals
 # plot(rep(randomeff[[1]][,1],each=length(residuals(a))/length(randomeff[[1]][,1])),
 #      residuals(a),xlab="random intercept",ylab="residual",
 #      main=paste(c("independence of random intercepts & \n residuals",inversion),collapse=" "))
 # # independence of random slopes & residuals
 # plot(rep(randomeff[[1]][,2],each=length(residuals(a))/length(randomeff[[1]][,2])),
 #      residuals(a),xlab="random slope",ylab="residual",
 #      main=paste(c("independence of random slopes & \n  residuals",inversion),collapse=" "))
 #  # check normal distribution of residuals
 #  qqnorm(residuals(a),main=paste(c("QQ Plot for residuals",inversion),collapse=" "))
 #  qqline(residuals(a),col="red")
  # check normal distribution random effects
 print(summary(a))
  qqnorm(ranef(a)[[1]][,1],main=paste(c("QQ Plot for random intercept locus",inversion),collapse=" "))
  qqline(ranef(a)[[1]][,1],col="red")
  qqnorm(ranef(a)[[2]][,1],main=paste(c("QQ Plot for random intercept replicates",inversion),collapse=" "))
  qqline(ranef(a)[[2]][,1],col="red")
  qqnorm(resid(a),main=paste(c("QQ Plot for residuals",inversion),collapse=" "))
  qqline(resid(a),col="red")
  plot(a)
}  
  
boxplot()

p.vals_v <- c()
p.vals_i <- c()

p.vals <- c()
for (inversion in names(a_trident) ) {
  a_nr <- lmer( arcsine(freq) ~ 1 + (1|repl)  +  pop +  pig, data = a_trident[[inversion]],control = lmerControl(optimizer ="Nelder_Mead"))
  #a1 <- lmer( arcsine(freq) ~ 0 + (1|pop/rep) + (0+pig|pop/rep) + pig, data = a_trident[[inversion]])
  #a <- lmer( arcsine(freq) ~ 1 + pop + (pig|repl) + pig, data = a_trident[[inversion]],control = lmerControl(optimizer ="Nelder_Mead"))
  a <- lmer( arcsine(freq) ~ 1 +  pop + (1 + pig |repl) + pig, data = a_trident[[inversion]],control = lmerControl(optimizer ="Nelder_Mead"))
  print(inversion)
  print(anova(a_nr,a,refit=FALSE))
  print(summary(a))
  p.vals=c(p.vals,anova(a)["pig","Pr(>F)"])
  names(p.vals)[length(p.vals)]=inversion
}
p.adjust(p.vals)
#> p.adjust(p.vals,"fdr")
#2L|t        2R|ns         3L|p         3R|p         3R|c        3R|mo 
#3.175860e-02 2.409389e-01 9.744444e-01 3.153071e-07 3.175860e-02 4.232152e-01 
p.vals <- c()
for (inversion in names(a_trident) ) {
  a <- lmer( arcsine(freq) ~ 1 + (1 + pig | repl) + pop  + pig, data = a_trident[[inversion]],control = lmerControl(optimizer ="Nelder_Mead"))
  p.vals=c(p.vals,anova(a)["pig","Pr(>F)"])
  names(p.vals)[length(p.vals)]=inversion
  randomeff <-ranef(a)
  # # independence of random intercept & random slope
  # plot(randomeff[[1]][,1],randomeff[[1]][,2],
  #      xlab ="random intercept" , ylab ="random slope",
  #      main=paste(c("correlated intercept & slope model",inversion),collapse=" "))
  # independence of random intercepts & residuals
  plot(rep(randomeff[[1]][,1],each=length(residuals(a))/length(randomeff[[1]][,1])),
       residuals(a),xlab="random intercept",ylab="residual",
       main=paste(c("independence of random intercepts & \n residuals",inversion),collapse=" "))
  # # independence of random slopes & residuals
  # plot(rep(randomeff[[1]][,2],each=length(residuals(a))/length(randomeff[[1]][,2])),
  #      residuals(a),xlab="random slope",ylab="residual",
  #      main=paste(c("independence of random slopes & \n  residuals",inversion),collapse=" "))
  # check normal distribution of residuals
  qqnorm(residuals(a),main=paste(c("QQ Plot for residuals",inversion),collapse=" "))
  qqline(residuals(a),col="red")
  # check normal distribution random effects
  qqnorm(ranef(a)[[1]][,1],main=paste(c("QQ Plot for random intercepts",inversion),collapse=" "))
  qqline(ranef(a)[[1]][,1],col="red")
  
}
p.adjust(p.vals,method="fdr")

  randomeff <-ranef(a)
  # # independence of random intercept & random slope
  # plot(randomeff[[1]][,1],randomeff[[1]][,2],
  #      xlab ="random intercept" , ylab ="random slope",
  #      main=paste(c("correlated intercept & slope model",inversion),collapse=" "))
  # independence of random intercepts & residuals
  plot(rep(randomeff[[1]][,1],each=length(residuals(a))/length(randomeff[[1]][,1])),
       residuals(a),xlab="random intercept",ylab="residual",
       main=paste(c("independence of random intercepts & \n residuals",inversion),collapse=" "))
  # # independence of random slopes & residuals
  # plot(rep(randomeff[[1]][,2],each=length(residuals(a))/length(randomeff[[1]][,2])),
  #      residuals(a),xlab="random slope",ylab="residual",
  #      main=paste(c("independence of random slopes & \n  residuals",inversion),collapse=" "))
  # check normal distribution of residuals
  qqnorm(residuals(a),main=paste(c("QQ Plot for residuals",inversion),collapse=" "))
  qqline(residuals(a),col="red")
  # check normal distribution random effects
  qqnorm(ranef(a)[[1]][,1],main=paste(c("QQ Plot for random intercepts",inversion),collapse=" "))
  qqline(ranef(a)[[1]][,1],col="red")
  
  #p.vals=c(p.vals,summary(a)[["coefficients"]][2,5])
  #names(p.vals)[length(p.vals)]=inversion
#  v <- lmer( arcsine(freq) ~ 1 + (1 |repl) + (0 + pig | repl) + pig, data = subset(a_trident[[inversion]], pop == "V"))
#  i <- lmer( arcsine(freq) ~ 1 + (1 |repl) + (0 + pig | repl) + pig, data = subset(a_trident[[inversion]], pop == "I"))
  p.vals=c(p.vals,anova(a)["pig","Pr(>F)"])
  names(p.vals)[length(p.vals)]=inversion
}
p.adjust(p.vals,method="bonferroni")


p.adjust(p.vals_v,method="bonferroni")
p.adjust(p.vals_i,method="bonferroni")



p.vals_v=c(p.vals_v,anova(v)[["Pr(>F)"]][1])
names(p.vals_v)[length(p.vals_v)]=inversion
p.vals_i=c(p.vals_i,anova(i)[["Pr(>F)"]][1])
names(p.vals_i)[length(p.vals_i)]=inversion
print(inversion)
print("Vienna")
print(anova(v))
print(inv.arcsine(summary(v)[["coefficients"]][,1]))
print("Italy")
print(anova(i))
print(inv.arcsine(summary(i)[["coefficients"]][,1]))

plot(a, resid(., scaled=TRUE) ~ fitted(.) | pig, abline = c(0,0))
qqmath(a)
plot(a)

boxplot( freq *100 ~ pig + repl   , data=a_trident[[inversion]],col=(c("yellow","brown")))
boxplot( arcsine(freq) ~ pig + repl   , data=a_trident[[inversion]],col=(c("yellow","brown")))


lmer_result_trident=list()
lmer_result_trident_vie=list()
lmer_result_trident_ita=list()
p.vals = c()
p.vals_v = c()
p.vals_i = c()
for (inversion in names(a_trident)){
  lmer_result_trident[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 
                                                1 + (1|pop/rep) + (1|locus)  + pig, data=a_trident[[inversion]], 
                                    family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  print(inversion)
  print(summary(lmer_result_trident[[inversion]])[["coefficients"]])
  p.vals=c(p.vals,summary(lmer_result_trident[[inversion]])[["coefficients"]][2,4])
  names(p.vals)[length(p.vals)]=inversion
  lmer_result_trident_vie[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 
                                                1 + (1|repl) + (1|locus) + pig, data=subset(a_trident[[inversion]], pop == "V"), 
                                              family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  p.vals_v=c(p.vals_v,summary(lmer_result_trident_vie[[inversion]])[["coefficients"]][2,4])
  names(p.vals_v)[length(p.vals_v)]=inversion
  
  lmer_result_trident_ita[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 
                                                1 + (1|repl) + (1|locus) + pig, data=subset(a_trident[[inversion]], pop == "I"), 
                                              family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  p.vals_i=c(p.vals_i,summary(lmer_result_trident_ita[[inversion]])[["coefficients"]][2,4])
  names(p.vals_i)[length(p.vals_v)]=inversion
  
}
inversion="3L|p"
inversion="3R|p"
inversion="2L|t"
inversion="3R|c"
inversion="3R|mo"
inversion="2R|ns"
for (inversion in names(a_trident)){
  boxplot( arcsine(freq) ~ pig + repl, data=a_trident[[inversion]],col=(c("yellow","brown")),main=inversion)
}
summary(glmer(cbind(inv_count, not_inv_count) ~ 
        1 + (1|pop) + (1|repl) + (1|locus) + pig, data=a_trident[[inversion]], 
      family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
summary(glmer(cbind(inv_count, not_inv_count) ~ 
                1 +   (1|repl) + (1|locus) + pig, data=subset(a_trident[[inversion]], pop == "I"), 
              family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
summary(glmer(cbind(inv_count, not_inv_count) ~ 
                1 + (1|repl) + (1|locus) + pig, data=subset(a_trident[[inversion]], pop == "V"), 
              family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
for (inversion in names(a_trident)){
  print(inversion)
  print(summary(lmer_result_trident_vie[[inversion]]))
}
for (inversion in names(a_trident)){
  print(inversion)
  print(summary(lmer_result_trident[[inversion]]))
}
for (inversion in names(a_trident)){
  print(inversion)
  print(summary(lmer_result_trident_ita[[inversion]]))
}

lmer_binom_result_11[[inversion]]<-lmer(cbind(inv_count, not_inv_count) ~ 
                                          1 + (1|prepl)+ (1|pop|repl/experiment) + pig, data=a_trident[[inversion]], family=binomial)

  
library(reshape2)

melt(a_11[[1]])

summaryBy(pop1 + pop2 + pop3 + pop4 +pop5 + pop6 ~ CHR + INV, data = inv_data_2010, FUN =c(median), keep.names = T)
#aggregate(cbind(pop1,pop2,pop3,pop4)~CHR+INV,data=inv_data,FUN=median)
a_10 <- splitBy(formula=~ CHR + INV, data = inv_data_2010 )
c_10 <- splitBy(formula=~ CHR + INV, data = inv_count_2010 )
for(i in names(c_10)){c_10[[i]] = sapply(c_10[[i]][,-(1:5)],median)}
for(i in c_10){ a = c(i[seq(from=1,to=length(i),by=2)]*i[seq(from=2,to=length(i),by=2)],(1-i[seq(from=1,to=length(i),by=2)])*i[seq(from=2,to=length(i),by=2)] ) }
for(i in c_10){ print(mantelhaen.test( array(matrix( c(i[seq(from=1,to=length(i),by=2)]*i[seq(from=2,to=length(i),by=2)],(1-i[seq(from=1,to=length(i),by=2)])*i[seq(from=2,to=length(i),by=2)]) ,nrow=4,byrow=T),dim=c(2,2,3)) )$p.value) }
c_10 <- splitBy(formula=~ CHR + INV, data = inv_count_2010 )
c_11 <- splitBy(formula=~ CHR + INV, data = inv_count_2011 )
for(i in names(c_11)){c_11[[i]] = sapply(c_11[[i]][,-(1:5)],median)}
for(i in c_11){ print(mantelhaen.test( array(matrix( c(i[seq(from=1,to=length(i),by=2)]*i[seq(from=2,to=length(i),by=2)],(1-i[seq(from=1,to=length(i),by=2)])*i[seq(from=2,to=length(i),by=2)]) ,nrow=4,byrow=T),dim=c(2,2,3)) )$p.value) }





quartz()
boxplot_inv(a_10,1,tit=", Vienna 2010")

quartz()
boxplot_inv(a_10,1,tit=", Vienna 2010")
dev.copy2pdf(file="inv_cs_Vie2010.pdf")
quartz()
boxplot_inv_combine(a_10,tit="Vienna 2010")
dev.copy2pdf(file="inv_cs_all_Vie2010.pdf")
quartz()
boxplot_inv_combine(a_11,tit="Bolzano 2011")
dev.copy2pdf(file="inv_cs_all_Ita2011.pdf")

inv_data_over = read.table("vie2010_2011.merged.inversions",header=T)
quartz()
boxplot_inv(a_10,6,tit=", Vienna 2010")
dev.copy2pdf(file="inv_3Rp_Vie2010.pdf")
for(i in 1:6){
	quartz()
boxplot_inv(a_11,i,tit=", Vienna 2010")
	}
# transforming list entries into matrices for marlies:
# count cov rep b_cs

fill_data <- function(raw_list){
	data_list = list()
	for(inversion in names(raw_list)) {
		int_tab=raw_list[[inversion]]
		for(i in seq(from=6,to=length(int_tab[1,]),by=2)){int_tab[,i]=int_tab[,i]*int_tab[,i+1]}
		marl_mat=matrix(nrow=0,ncol=4)
	colnames(marl_mat)=c("count","cov","rep","cs")
	for(i in 1:3){
		j=(i-1)*2
		loc_mat=cbind(int_tab[,(6+j):(6+j+1)],i,1)
		colnames(loc_mat)=c("count","cov","rep","cs")
		marl_mat=rbind(marl_mat,loc_mat) 
		loc_mat=cbind(int_tab[,(12+j):(12+j+1)],i,0)
		colnames(loc_mat)=c("count","cov","rep","cs")
		marl_mat=rbind(marl_mat,loc_mat)
		}
	rownames(marl_mat)=NULL
	data=as.data.frame(marl_mat)	
	data$freq<-round((data$count/data$cov)*100,digits=0)
	data$inv_count<-round(data$count,digits=0)
	data$not_inv_count<- data$cov-data$inv_count
	data$experiment<-as.factor(data$rep*10+data$cs)
	data$rep_fac<-as.factor(data$rep)
	data$cs_fac<-as.factor(data$cs)
	data_list[[inversion]] = data		
	}		
	return(data_list)	
}	

data_11=fill_data(a_11)	
data_10=fill_data(c_10)	
str(data_10)
library("lme4")
lmer_result_11=list()
for (inversion in names(data_11)){
	lmer_result_11[[inversion]]<-lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_11[[inversion]], family=poisson)
}
lmer_result_10=list()
for (inversion in names(data_10)){
	lmer_result_10[[inversion]]<-lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10[[inversion]], family=poisson)
}
summary(lmer_result_10[[1]])
lmer_binom_result_11=list()
for (inversion in names(data_11)){
	lmer_binom_result_11[[inversion]]<-lmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_11[[inversion]], family=binomial)
}
summary(lmer_binom_result_11[[1]])
lmer_binom_result_10=list()
for (inversion in names(data_10)){
	lmer_binom_result_10[[inversion]]<-lmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10[[inversion]], family=binomial)
}
summary(lmer_binom_result_10[[1]])

# only significant result in Bolzano, no significant result for in(2L)t in Vienna
# try 10 2L_t without replicate 2
data_10_2Lt_wo_rep2=data_10[["2L|t"]][data_10[["2L|t"]]$rep != 2, ]
summary(lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10_2Lt_wo_rep2, family=poisson))
summary(lmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10_2Lt_wo_rep2, family=binomial))
# install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R", type="source")
library(coefplot2)
coefplot2(test_lmer)
coeftab(test_lmer)
library(influence.ME)
x11()
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
f <- fitted(fm1)
r <- residuals(fm1)
plot(f,r)
sm <- loess(r~f)
v <- seq(min(f),max(f),length=101)
lines(v,predict(sm,data.frame(f=v)),col=2)
test_lmer=lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10_2Lt_wo_rep2,family=poisson)
f <- fitted(test_lmer)
r <- residuals(test_lmer)
plot(f,r)
sm <- loess(r~f)
v <- seq(min(f),max(f),length=101)
lines(v,predict(sm,data.frame(f=v)),col=2)
library(ggplot2)
ggplot(data = sleepstudy, aes(x = Days, y = Reaction, color = factor(Subject))) + geom_line(aes(group = Subject)) + geom_point()
