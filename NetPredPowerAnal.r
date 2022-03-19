## R code for PO1 outcome analysis ##

## Longitudinal analysis yields more interesting results; revisit appears to matter! ##

#Input data sets contain t-statistics of network predictors at predicting ICD status and overall c-stat
#Separate coefficients and c-statistic for each design combination
#Invert to find design that yields coefficient of largest magnitude for each network statistic
# Models with directed measures: AdoptDirMeasures2011lnrAll_clean.csv
# Models with undirected measures: AdoptUndirMeasures2011lnrAll_clean.csv

#Data Dictionary (just key variables used so far) for hrr_network_stats.csv:
#Interm = 0 implies that ABC yields edges AB and AC only (Interm = 1 adds AC: allow interm encounter)?
# Note: I now reverse the sign of Interm to conform with Continuity
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
library(gee)
library(gam)

#Set working directory (you need to set this up); I use a server (substitute your directory)
rsource='/Volumes/rc/Dartmouth/Biostatistics/P01/Code'
setwd(rsource) #Working directory
datdir='../Data/' #Directory to store data (you need to set this up)
outdir='../Output/' #Directory to store output (you need to set this up)

#Load ICC summary measures: Directed
indata=paste(datdir,'AdoptDirMeasures2011lnrAll_clean.csv',sep="")
dataicd=data.frame(read.csv(file=indata,header=T))
names(dataicd)
dataicd$Interm=1-dataicd$Interm #So that continuity enforced when = 1

#Model of c-statistic regressed on network design factors
modoptdir=summary(lm(cprob~factor(Type)+Interm+Revisit+Binary,data=dataicd))$coef

#Just keep Type 1 design for binary-valued designs to emulate icc/variation-based analyses
trimweight1 <- 0
if (trimweight1==1) {
  ind <- (dataicd$Type==1)
  dataicd=dataicd[ind,]
}

#Top 10 measures based on decrease from baseline in deviance, stratified by whether network is binary as diff base used
ind <- (dataicd$Binary==0)
dataicdw <- dataicd[ind,]
ord <- order(dataicdw$cprobDiff,decreasing=TRUE)
top10wdir <- data.frame(dataicdw[ord[1:5],c(1:6,23:30)]) #Features of best design
dataicdb <- dataicd[!ind,]
ord <- order(dataicdb$cprobDiff,decreasing=TRUE)
top10bdir <- data.frame(dataicdb[ord[1:5],c(1:6,23:30)]) #Features of best design
top10dir <- rbind(top10wdir,top10bdir)

#Best construction for making each network measure the most predictive
nmeasure=19 #columns 7 thru 23, 25, 27, 29, 30
cols <- c(7:23,25,27,29,30)
bestdsgn <- apply(dataicd[,cols],2,function(x) order(abs(x),decreasing=TRUE)[1])
maxz <- apply(dataicd[,cols],2,function(x) max(abs(x)))
maxzsign <- diag(as.matrix(dataicd[bestdsgn,cols]))
optdsgn10dir <- data.frame(dataicd[bestdsgn,1:6],maxz,round(maxzsign,3)) #Features of best design
row.names(optdsgn10dir) <- names(maxz)

#Optimal fitted models of weighted and binary directed designs: Used in Table 6
ind <- (dataicd$Mthd %in% c(20,42))
modoptdir <- dataicd[ind,]


#Load ICC summary measures: Undirected is meaningless so don't use
indata=paste(datdir,'AdoptUndirMeasures2011lnrAll_clean.csv',sep="")
dataicdu=data.frame(read.csv(file=indata,header=T))
names(dataicdu)
dataicdu$Interm=1-dataicdu$Interm

#Model of c-statistic regressed on network design factors
modoptun=summary(lm(cprob~factor(Type)+Interm+Revisit+Binary,data=dataicdu))$coef

trimweight1 <- 0
if (trimweight1==1) {
  ind <- (dataicdu$Type==1)
  dataicdu=dataicdu[ind,]
}

#Top 10 measures based on c-statistic
ind <- (dataicdu$Binary==0)
dataicduw <- dataicdu[ind,]
ord <- order(dataicduw$cprobDiff,decreasing=TRUE)
top10wudir <- data.frame(dataicduw[ord[1:5],c(1:6,15:22)]) #Features of best design
dataicdub <- dataicdu[!ind,]
ord <- order(dataicdub$cprobDiff,decreasing=TRUE)
top10budir <- data.frame(dataicdub[ord[1:5],c(1:6,15:22)]) #Features of best design
top10udir <- rbind(top10wudir,top10budir)

#Best construction for each network measure the most predictive
cols <- c(7:15,17,19,21,22)
bestdsgn <- apply(dataicdu[,cols],2,function(x) order(abs(x),decreasing=TRUE)[1])
maxz <- apply(dataicdu[,cols],2,function(x) max(abs(x)))
maxzsign <- diag(as.matrix(dataicdu[bestdsgn,cols]))
optdsgn10udir <- data.frame(dataicdu[bestdsgn,1:6],maxz,maxzsign) #Features of best design
row.names(optdsgn10udir) <- names(maxz)

optudirdsgn <- c(5,23) #Weighted and Binary optimal design based on change in c-statistic
modoptudir <- dataicdu[optudirdsgn,]

#Base fitted models (26 for weighted, 66 for binary)
ind <- (dataicdu$Mthd %in% c(26,66))
modbase <- dataicdu[ind,]


#Save output used in Tables 5 and 6
sink(paste(outdir,'ICDOutcomeResults2011lnrAll.csv',sep=""))
print("Top ranked directed designs")
print(round(top10dir,digits=4))
print("Top ranked directed designs by feature")
print(round(optdsgn10dir[,c(2:4,6,8)],digits=4))
print("Coefficients directed designs")
print(round(modoptdir[,c(2:22)],digits=4))

print("Top ranked undirected designs")
print(round(top10udir,digits=4))
print("Top ranked undirected designs by feature")
print(round(optdsgn10udir[,c(2:4,6,8)],digits=4))
print("Coefficients undirected designs")
print(round(modoptudir[,c(2:14)],digits=4))

print("Printed models for base (undirected) designs")
print(round(modbase[,c(2:14)],digits=4))
sink()

