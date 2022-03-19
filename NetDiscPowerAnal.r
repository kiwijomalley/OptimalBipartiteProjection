## R code for PO1 exploratory analysis ##

#Data sets:
# ICCnetPHNall2010nr.csv = Data set with intraclass correlation coefficients
#  = Data set with network-level summary measures

#Outputs:
# ICCResults2010nr.csv
# Contains Results appearing in Tables 2 through 4 and additional output

#Data Dictionary (just key variables used so far) for hrr_network_stats.csv:
#Interm = 0 implies that ABC yields edges AB and AC only (Interm = 1 adds AC: allow interm encounter)?
#Revisit = 1 implies revisit constraint is required (ABA must occur for contribution to AB edge)? 
#Directed = 1 implies directed network
#Binary = 1 implies binary network
#Weights: 1 = W1, 2 = W2, 3 = W3, 4 = W2 (normalized), 5 = W3 (normalized)
#Most restricted combination is Interm = 0 and Revisit = 1 (must revisit with no interim encounter)

#Load libraries
library(reshape)
library(date)
library(Hmisc)
library(lme4)
library(Zelig)
#library(ZeligMultilevel)
library(gee)
library(gam)

#Set working directory (you need to set this up); I use a server (substitute your directory)
rsource='/Volumes/OMalleyJ/Dartmouth/Biostatistics/P01/Code'
setwd(rsource) #Working directory
datdir='../Data/' #Directory to store data (you need to set this up)
outdir='../Output/' #Directory to store output (you need to set this up)

#Load ICC summary measures
indata=paste(datdir,'ICCnetPHNall2010nr.csv',sep="")
#indata=paste(datdir,'ICCnetPHNall2010sc.csv',sep="") #Old
dataicc=data.frame(read.csv(file=indata,header=T))
names(dataicc)
dataicc$Type=rep(seq(1,5),times=80/5) #Type of weights (W1, W2, W3, W2n, W3n)
dataicc$Interm=rep(rep(c(0,1),each=5),times=80/10) #1 = Intermediate encounters allowed, 0 = Not
dataicc$Revisit=rep(rep(c(0,1),each=10),times=80/20) #1 = Revisit required, 0 = Not
dataicc$Directed=rep(rep(c(1,0),each=20),times=80/40) #1 = Directed, 0 = Undirected
dataicc$Binary=rep(rep(c(0,1),each=40),times=80/80) #1 = Binary, 0 = Non-binary
dataicc=dataicc[,c(1,7:11,2:6)]
dataicc$Interm=1-dataicc$Interm

#Load network level data for each measure
indata=paste(datdir,'ICCnetmeaPHN2010nr.csv',sep="")
#indata=paste(datdir,'ICCnetmeaPHN2010sc.csv',sep="")
data=data.frame(read.csv(file=indata,header=T))
names(data)
nmeasure=16
data$Type=rep(rep(seq(1,5),times=80/5),each=nmeasure)
data$Interm=rep(rep(rep(c(0,1),each=5),times=80/10),each=nmeasure)
data$Revisit=rep(rep(rep(c(0,1),each=10),times=80/20),each=nmeasure)
data$Directed=rep(rep(rep(c(1,0),each=20),times=80/40),each=nmeasure)
data$Binary=rep(rep(rep(c(0,1),each=40),times=80/80),each=nmeasure)
data=data.frame(data)
MeasNames=c("Nodes","Isolates","Density","Specialization","CorDeg_InOut","Assort_InIn","Assort_InOut","Assort_OutIn","Assort_OutOut","Reciprocity","Transitivity","Ave_Clustering","Ncomponents","Diameter","sumTransitivity","sum3Cycle")
data$MeasureName=rep(MeasNames,times=80) #sum terms add up over triad census
data=data[,c(1,10,2,5:9,4,3)]
#Drop measures that are not self-scaled (e.g., proportions) or invariant to scale
ind <- (data$Measure==3) #Check if density varies so much between designs that was unscaled
tmp <- data[ind,]
tapply(tmp$wvar,tmp$Type,'mean',na.rm=TRUE)
dropy = 0
if (dropy==1) { 
 dropMeas=c(3,4)
 ind <- (data$Measure %in% dropMeas)
 data <- data[!ind,]
}
data$Interm=1-data$Interm #So that 

sink(paste(outdir,'ICCResults2010nr.csv',sep=""))

print("Form Unique Designs of Network Summary Statistics")
#Only 3 factors should vary when network is binary if have 0 threshold 
# These factors are: Intermediate, Revisit, Directed
ind=(dataicc$Binary==1)
dataiccb=dataicc[ind,]
ind=(dataicc$Binary==0)
dataiccw=dataicc[ind,]

#Analyzed data sets (a=all, u=undirected, d=directed)
dataicca=rbind(dataiccw,dataiccb) #Very little variation among weighting schemes
ind=(dataicca$Directed==0)
dataiccu=dataicca[ind,]
dataiccd=dataicca[!ind,]

#Rank the methods of building the network overall based on icc (also include the worst measure)
# Used in network-wide part of Table 2 and Table 3
print("Analyze undirected measures over undirected measures (same as over all measures)")
ord=order(dataiccu$iccu,decreasing = TRUE)
print("Top 10 undirected designs based on overall icc model")
print(dataiccu[ord[1:10],])
modiccu=lm(iccu~factor(Type)+Interm+Revisit+Binary,data=dataiccu)
print("Main Effects Model")
print(summary(modiccu))

print("Analyze directed measures over all measures")
ord=order(dataiccd$icc,decreasing = TRUE)
print("Top 10 directed designs based on overall icc model")
print(dataiccd[ord[1:10],])
modiccd=lm(icc~factor(Type)+Interm+Revisit+Binary,data=dataiccd)
print("Key output: Main Effects Model")
print(summary(modiccd))

#Extension to consider interaction effects
print("Analyze directed measures over all measures with interaction effects for the primary factors")
ord=order(dataiccd$icc,decreasing = TRUE)
print("Top 10 directed designs based on overall icc model")
print(dataiccd[ord[1:10],])
modiccdint=lm(icc~factor(Type)+Interm+Revisit+Binary+Interm*Revisit+Interm*Binary+Revisit*Binary+Interm*Revisit*Binary,data=dataiccd)
print("Key output: Interaction Effects Model")
print(summary(modiccdint))
modiccdintr=lm(icc~factor(Type)+Interm+Revisit+Binary+Interm*Revisit,data=dataiccd)
print("Key output: Reduced Interaction Effects Model")
print(summary(modiccdintr))
print(anova(modiccd,modiccdintr,modiccdint))

print("Analyze directed measures over directed measures")
ord=order(dataiccd$iccd,decreasing = TRUE)
print("Top 10 directed designs based on overall icc model for just the directed measures")
print(dataiccd[ord[1:10],])
modiccd=lm(iccd~factor(Type)+Interm+Revisit+Binary,data=dataiccd)
print("Main Effects Model")
print(summary(modiccd))


## Measure-specific measures: Used in Table 3 ##

print("Form Unique Designs of Network Summary Statistics")
#Only 3 factors should vary when network is binary if have 0 threshold 
# These factors are: Intermediate, Revisit, Directed
#ind=(data$Binary==1 & data$Type==1) #Setting Type = 1 extracts first occurrence of binary design as weights immaterial
ind=(data$Binary==1)
datab=data[ind,]
ind=(data$Binary==0)
dataw=data[ind,]

#Analyzed data sets (a=all, u=undirected, d=directed)
dataa=rbind(dataw,datab) #Very little variation among weighting schemes
ind=(dataa$Directed==0)
datau=dataa[ind,]
datad=dataa[!ind,]
#Because weighting had no impact on undirected network, remove all cases for which Type is other than 1
ind=(datau$Type==1)
datau=datau[ind,]

print("Network Discrimination of Undirected Measures over undirected Network designs")
UnDir=c(1:4,11:14)
if (dropy==1) {
 UnDir=c(1,2,11:14)
}
ind=(datau$Measure %in% UnDir)
datauu=datau[ind,]
nmthd <- length(unique(datauu$Mthd))
nu=length(UnDir)
val=matrix(0,nrow=nmthd,ncol=nu)
rankmthd=matrix(0,nrow=nmthd,ncol=nu)
bestdsgnu <- data.frame(matrix(0,nrow=nu,ncol=ncol(datauu)))
fitmthd=matrix(0,nrow=nu,ncol=4)
rankall=matrix(0,nrow=nmthd,ncol=nu)
for (i in 1:nu) {
 ind=(datauu$Measure==UnDir[i])
 tmp=datauu[ind,]
 val[,i]=tmp$wvar
 rankmthd[,i]=rank(-tmp$wvar,ties.method="average") #Gives ranking of method on each given network measure
 pos=order(tmp$wvar,decreasing=TRUE)[1]
 bestdsgnu[i,] <- tmp[pos,]
 mod=lm(wvar~Interm+Revisit+Binary,data=tmp)
 fitmthd[i,]=summary(mod)$coef[,3]
}
rankave <- apply(rankmthd,1,'mean')
rankdsgn <- data.frame(cbind(unique(datauu$Mthd),rankmthd,rankave))
names(rankdsgn) <- c("Mthd",MeasNames[UnDir],"AveRank")
val=data.frame(val)
valmn <- apply(val,2,'mean')
augvalmn <- rep(valmn,times=nmthd)
datauu$swvar <- datauu$wvar/augvalmn
#tapply(datauu$swvar,datauu$Measure,'mean')
rankmat <- as.matrix(rankdsgn[,c(2:(nu+1))])
rankvec <- data.frame(values=c(t(data.matrix(rankmat))), ind=I(rep(colnames(rankmat), nrow(rankmat))))
datauu$srank <- rankvec$values
#tapply(datauu$srank,datauu$Measure,'mean')

print("Undirected Measures")
valstd <- val / (as.matrix(rep(1,nrow(val)),ncol=1) %*% t(as.matrix(valmn)))
rankdsgn$stdaveval <- apply(valstd,1,'mean')
#rankdsgn <- data.frame(dataiccu[(dataiccu$Type==1),c(1:4,6)], rankdsgn[,-1])
rankdsgn <- data.frame(dataiccu[,c(1:4,6)], rankdsgn[,-1])
ord <- order(rankdsgn$AveRank)
rankdsgn <- rankdsgn[ord,]
print("Designs ordered by average rank with measure-specific results and also average standardized within-variance")
print(rankdsgn)
names(bestdsgnu) <- names(datauu)[1:(ncol(datauu)-2)]
print("Best design for each measure")
print(bestdsgnu)
fitmthd=data.frame(Measure=MeasNames[UnDir],fitmthd)
names(fitmthd)[2:5]=names(mod$coef)
print("Fits for undirected designs across undirected measures")
print(fitmthd)
#Final regression model that pools across measures after standardizing wvar by measure
datauu$lswvar<-log(datauu$swvar+0.05)
modsu=lm(swvar~Interm+Revisit+Binary,data=datauu)
print("Overall regression model for undirected case")
print(summary(modsu))
#Final regression model that pools across measures and regresses ranks on design factors
modru=lm(srank~Interm+Revisit+Binary,data=datauu)
print("Overall regression model of srank for undirected case")
print(summary(modru))

print("Network Discrimination of All Measures in Directed Network: Used in Table 4")
Dir=1:16 #Want over all measures as all measures are well-defined when network is directed
if (dropy==1) {
 Dir=c(1,2,5:16)
}
#Dir=c(5:10,15:16) 
#ind=(data$Measure %in% Dir & data$Directed==1)
## Could consider restricting data to Type==1 for undirected measures (as their values don't seem to vary with Type)
datadd=datad
nmthd <- length(unique(datadd$Mthd))
nu=length(Dir)
val=matrix(0,nrow=nmthd,ncol=nu)
rankmthd=matrix(0,nrow=nmthd,ncol=nu)
bestdsgnd <- data.frame(matrix(0,nrow=nu,ncol=ncol(datadd)))
bestdsgnd2 <- data.frame(matrix(0,nrow=nu,ncol=ncol(datadd))) #Second best design
fitmthd=matrix(0,nrow=nu,ncol=8)
fitmthdint=matrix(0,nrow=nu,ncol=9)
rankall=matrix(0,nrow=nmthd,ncol=nu)
for (i in 1:nu) {
 ind=(datadd$Measure==Dir[i])
 tmp=datadd[ind,]
 val[,i]=tmp$wvar
 rankmthd[,i]=rank(-tmp$wvar,ties.method="average") #Gives ranking of method on each given network measure with lower rank = better
 pos2=order(tmp$wvar,decreasing=TRUE)[1:2] #Best two designs
 pos <- pos2[1] #Best design
 bestdsgnd[i,] <- tmp[pos,]
 bestdsgnd2[i,] <- tmp[pos2[2],]
 mod=lm(wvar~factor(Type)+Interm+Revisit+Binary,data=tmp)
 fitmthd[i,]=summary(mod)$coef[,3]
 modint=lm(wvar~factor(Type)+Interm+Revisit+Binary+Interm*Revisit,data=tmp) #Over-fits for many measures!
 fitmthdint[i,]=summary(modint)$coef[,3]
}
rankave <- apply(rankmthd,1,'mean')
rankdsgn <- data.frame(cbind(unique(datadd$Mthd),rankmthd,rankave))
names(rankdsgn) <- c("Mthd",MeasNames[Dir],"AveRank")
val=data.frame(val)
valmn <- apply(val,2,'mean')
augvalmn <- rep(valmn,times=nmthd)
datadd$swvar <- datadd$wvar/augvalmn
#tapply(datadd$swvar,datadd$Measure,'mean')
rankmat <- as.matrix(rankdsgn[,c(2:(nu+1))])
rankvec <- data.frame(values=c(t(data.matrix(rankmat))), ind=I(rep(colnames(rankmat), nrow(rankmat))))
datadd$srank <- rankvec$values
#tapply(datadd$srank,datadd$Measure,'mean')

print("Directed Measures") #Used in summary measure part of Table 4
valstd <- val / (as.matrix(rep(1,nrow(val)),ncol=1) %*% t(as.matrix(valmn)))
rankdsgn$stdaveval <- apply(valstd,1,'mean')
rankdsgn <- data.frame(dataiccd[,c(1:4,6)], rankdsgn[,-1])
ord <- order(rankdsgn$AveRank)
rankdsgn <- rankdsgn[ord,]
print("Designs ordered by average rank with measure-specific results and also average standardized within-variance")
print(rankdsgn)
names(bestdsgnd) <- names(datadd)[1:(ncol(datadd)-2)]
names(bestdsgnd2) <- names(datadd)[1:(ncol(datadd)-2)]
bestdsgnd$secvar <- round(100*(bestdsgnd$wvar - bestdsgnd2$wvar)/bestdsgnd$wvar,2)
print("Best design for each measure")
print(bestdsgnd[,c(1:9,11)])
print('Best design based on pooled standardized variance')
poolswvar <- tapply(datadd$swvar,datadd$Mthd,'mean')
pos <- order(poolswvar,decreasing=TRUE)[1:2]
print(poolswvar[pos]) 
ind <- (datadd$Mthd %in% c(13,18)) #Need to manually extract Mthd corresponding to maximum value
optswvardsgn <- datadd[ind,]
print(optswvardsgn[c(1,17),c(1,3:8)]);
print('Best design based on pooled standardized rank')
poolsrank <- tapply(datadd$srank,datadd$Mthd,'mean')
pos <- order(poolsrank,decreasing=FALSE)[1:2]
print(poolsrank[pos])
ind <- (datadd$Mthd %in% c(13,12)) #Need to manually extract Mthd corresponding to minimum value (lowest rank)
optsrankdsgn <- datadd[ind,]
print(optsrankdsgn[c(1,17),c(1,3:8)])
print("Fits for directed designs across all measures")
fitmthd=data.frame(Measure=MeasNames[Dir],fitmthd)
names(fitmthd)[2:9]=names(mod$coef)
print(fitmthd)
#Final regression model that pools across measures and regresses standardizing wvar on design factors
datadd$lswvar<-log(datadd$swvar+0.05)
modsd=lm(swvar~factor(Type)+Interm+Revisit+Binary,data=datadd)
print("Key output: Overall regression model of swvar for directed case")
print(summary(modsd))
#Final regression model that pools across measures and regresses ranks on design factors
modrd=lm(srank~factor(Type)+Interm+Revisit+Binary,data=datadd)
print("Key output: Overall regression model of srank for directed case")
print(summary(modrd)) #Coefficients need to be reverse!

#Interaction effect analyses of Final regression model that pools across measures and regresses standardizing wvar on design factors
datadd$lswvar<-log(datadd$swvar+0.05)
modsdint=lm(swvar~factor(Type)+Interm+Revisit+Binary+Interm*Revisit+Interm*Binary+Revisit*Binary+Interm*Revisit*Binary,data=datadd)
print("Key output: Overall regression model of swvar for directed case with interactions")
print(summary(modsdint))
modsdintr=lm(swvar~factor(Type)+Interm+Revisit+Binary+Interm*Revisit,data=datadd)
print(summary(modsdintr))
print(anova(modsd,modsdintr,modsdint))
#Final regression model that pools across measures and regresses ranks on design factors
modrdint=lm(srank~factor(Type)+Interm+Revisit+Binary+Interm*Revisit+Interm*Binary+Revisit*Binary+Interm*Revisit*Binary,data=datadd)
print("Key output: Overall regression model of srank for directed case with interactions")
print(summary(modrdint)) #Coefficients need to be reverse!
modrdintr=lm(srank~factor(Type)+Interm+Revisit+Binary+Interm*Revisit,data=datadd)
print(summary(modrdintr)) #Coefficients need to be reverse!
print(anova(modrd,modrdintr,modrdint))

#Expect negative correlation of these:
cor(cbind(datauu$ave,datauu$wvar,datauu$swvar,datauu$lswvar,datauu$srank))
cor(cbind(datadd$ave,datadd$wvar,datadd$swvar,datadd$lswvar,datadd$srank))

sink()