#Evaluation of predictive accuracy for predicting the likelihood of the PHN
# having adopted ICD capability. 
#We compare the level of predictability across the 80 designs (separate models per design)
#Only have phn data in 2011 so regress 2011 icd status on 2010 network measures

library(lme4)
setwd("/projects/active/14593/programs/jomalley/SmartNetwork")
datdir="/projects/active/14593/idata/core_c/core_c_jomalley/PHN-subnetwork/2010/output-"
phndir="/projects/active/14593/idata/core_c/core_c_jomalley/PHN-subnetwork/2011"
outdir="/projects/active/14593/idata/core_c/core_c_jomalley/PHN-subnetwork/2011/SmartNetwork/"

##########################################################################################
# Function OptimisedConc : for concordance, discordance, ties
# The function returns Concordance, discordance, and ties
# by taking a glm binomial model result as input.
# Although it still uses two-for loops, it optimises the code
# by creating initial zero matrices
# Reference: http://shashiasrblog.blogspot.com/2014/01/binary-logistic-regression-on-r.html
###########################################################################################
OptimisedConc=function(model)
{
  Data = cbind(model$y, model$fitted.values) 
  ones = Data[Data[,1] == 1,]
  zeros = Data[Data[,1] == 0,]
  conc=matrix(0, dim(zeros)[1], dim(ones)[1])
  disc=matrix(0, dim(zeros)[1], dim(ones)[1])
  ties=matrix(0, dim(zeros)[1], dim(ones)[1])
  for (j in 1:dim(zeros)[1])
  {
    for (i in 1:dim(ones)[1])
    {
      if (ones[i,2]>zeros[j,2])
      {conc[j,i]=1}
      else if (ones[i,2]<zeros[j,2])
      {disc[j,i]=1}
      else if (ones[i,2]==zeros[j,2])
      {ties[j,i]=1}
    }
  }
  Pairs=dim(zeros)[1]*dim(ones)[1]
  PercentConcordance=(sum(conc)/Pairs)*100
  PercentDiscordance=(sum(disc)/Pairs)*100
  PercentTied=(sum(ties)/Pairs)*100
  return(list("Percent Concordance"=PercentConcordance,"Percent Discordance"=PercentDiscordance,"Percent Tied"=PercentTied,"Pairs"=Pairs))
}

Nbatch=30 #30 when all data is produced
Nnet=160
posempty=c()
for (j in 1:Nbatch) {
 indir=paste(datdir,j,sep="")
 pos1=(j-1)*Nnet
 for (i in pos1:(pos1+Nnet-1)) {
  infile=paste(paste("PHN-",i,sep=""),"-network-measures.txt",sep="") #Chuankai
  f=try(read.table(paste(indir,infile,sep="/"),header=FALSE,sep=";"))
  if (class(f)!='try-error') {
   data <- read.table(paste(indir,infile,sep="/"),header=FALSE,sep=",")
   data=data[,1:(ncol(data)-1)]
   data$netid=i
   if (i==0) {
    netdata=data
   } else {
    netdata=rbind(netdata,data) 
   }
  } else {
   posempty=c(posempty,i)
  }
 }
}
Nnetwork=Nbatch*Nnet

names(netdata) <- c("size","nisolate","density_str","centralization","cor_inout","assort_inin","assort_inout","assort_outin","assort_outout","recip","trans","avg_clust","tc2","tc3","tc4","tc5","tc6","tc7","tc8","tc9","tc10","tc11","tc12","tc13","tc14","tc15","tc16","Ncomponent","diam","sumTrans","sum3Cycle","netid") #Temporary use until 2nd column figured out!

#Note: Size = non-isolates

#Set missing values
ind=(netdata==-99)
netdata[ind]=NA
#Due to Chuankai's code not computing measures when network has size<=20 modes all these network measures
# equal 0. Want to set the measures other than size and number of nodes to NA
#Consider getting this redone (maybe get errors when compute network measures)
indrow=(netdata$size<=20) 
netdata[indrow,3:(ncol(netdata)-1)]=NA #Don't want to override netid

#Add design variables for each of the Nnetworks
type=rep(seq(1,5),times=80/5)
interm=rep(rep(c(0,1),each=5),times=80/10)
revisit=rep(rep(c(0,1),each=10),times=80/20)
directed=rep(rep(c(1,0),each=20),times=80/40)
binary=rep(rep(c(0,1),each=40),times=80/80)
netdata$Mthd=rep(seq(1,80),times=Nnetwork) #Because start count at 1 don't need to use sMthd==(i-1) later!
netdata$Type=rep(type,times=Nnetwork)
netdata$Interm=rep(interm,times=Nnetwork)
netdata$Revisit=rep(revisit,times=Nnetwork)
netdata$Directed=rep(directed,times=Nnetwork)
netdata$Binary=rep(binary,times=Nnetwork)

wantstand <- 1 #Standardizing seems to lower c-stat
if (wantstand==1) {
 #Statistical transformations that help with comparability of the same measure across the 80 designs
 totnodes <- netdata$size+netdata$nisolate
 tottriads <- totnodes*(totnodes-1)*(totnodes-2)/6
 alltriads <- apply(netdata[,13:27],1,'sum') #Sum's to 100000 as samples were taken as opposed to evaluating over the entire network
 netdata$sumTrans <- netdata$sumTrans/100000 #Not perfect; get some values > 1
 netdata$sum3Cycle <- netdata$sum3Cycle/100000 #Not perfect; get some values > 1
 netdata$centralization <- netdata$centralization/((netdata$density_str*(netdata$size-1))^2)
 
relwtsbase <- 66 #Base = binary, undirected, type 1, intermediate allowed and revisit not required
 relwts <- tapply(netdata$density_str,netdata$Mthd,'mean',na.rm=TRUE)
 relwtsvec <- rep(relwts/relwts[relwtsbase],times=nrow(netdata)/80) #80 = Number of distinct designs
 netdata$density_str=netdata$density_str/relwtsvec #Makes density/strength have same mean as base design 
 tapply(netdata$density_str,netdata$Mthd,'mean',na.rm=TRUE) #Check
 tapply(netdata$density_str,netdata$Type,'mean',na.rm=TRUE) #Compare

 #Standardize network measures to make them comparable across measures (for pooled across measures analysis) by having SD of 1 across entire data set
 stdev <- sqrt(apply(netdata[,1:31],2,'var',na.rm=TRUE))
 ones <- rep(1,nrow(netdata))
 scalemat <- matrix(ones,ncol=1) %*% t(matrix(stdev,ncol=1))
 netdata[,1:31] <- netdata[,1:31]/scalemat
}

basemeasure <- c("size","nisolate","density_str","centralization","trans","avg_clust","Ncomponent","diam","netid","Binary")
indbase <- (netdata$Type==1 & netdata$Interm==1 & netdata$Revisit==0 & netdata$Directed==0) #Mthd = 26 (weighted), 66 (binary)
datbase <- netdata[indbase,basemeasure]
names(datbase)[1:8] <- paste(basemeasure[1:8],"b",sep="")
netdata <- merge(netdata,datbase,by=c("netid","Binary"))

#Link network data to phn icd data via crosswalk of phn IDs
linkphn <- read.table(paste(phndir,"phn_CAind_comp_2011_v2.csv",sep="/"),header=TRUE,sep=",") #Only have 2011 data
icdphn <- read.table(paste(phndir,"ICDegoalter.csv",sep="/"),header=TRUE,sep=",")
linkout <- merge(linkphn,icdphn[(icdphn$year==5),],by="phn")
analdata <- merge(netdata,linkout,by.x="netid",by.y="phn_index")

#Trim data to just hospitals for which they have non-missing data across all designs (get from Directed=1, Revisit=1, and Interm=0)
missdsgn <- (analdata$Directed==1 & analdata$Revisit==1 & analdata$Interm==0 & analdata$Binary==0 & analdata$Type==1) #The design most prone to missing values
indmissdir <- is.na(analdata$sizeb+analdata$nisolateb+analdata$Ncomponentb+analdata$diamb+analdata$density_strb+analdata$centralizationb+analdata$transb+analdata$avg_clustb+analdata$size+analdata$nisolate+analdata$Ncomponent+analdata$diam+analdata$density_str+analdata$centralization+analdata$cor_inout+analdata$assort_inin+analdata$assort_inout+analdata$assort_outin+analdata$assort_outout+analdata$recip+analdata$trans+analdata$avg_clust+analdata$sumTrans+analdata$sum3Cycle)
goodhosp <- unique(analdata$netid[(missdsgn & !indmissdir)])
ind <- (analdata$netid %in% goodhosp)
analdata <- analdata[ind,] #Ensures that use same set of hospitals for all analyses

#Evaluate which undirected measures vary from their baseline value (Result: Nothing is invariant)
sizecheck <- analdata$size-analdata$sizeb; summary(sizecheck)
pisocheck <- analdata$nisolate-analdata$nisolateb; summary(pisocheck)
tsize <- analdata$size+analdata$nisolate;
tsizeb <- analdata$sizeb+analdata$nisolateb;
tsizecheck <- tsize-tsizeb; summary(tsizecheck)
Ncomcheck <- analdata$Ncomponent-analdata$Ncomponentb; summary(Ncomcheck)
diamcheck <- analdata$diam-analdata$diamb; summary(diamcheck)
dencheck <- analdata$density_str-analdata$density_strb; summary(dencheck)
cencheck <- analdata$centralization-analdata$centralizationb; summary(cencheck)
trancheck <- analdata$trans-analdata$transb; summary(trancheck)
cluscheck <- analdata$avg_clust-analdata$avg_clustb; summary(cluscheck)

#Models of ICD status (stratified by method)
nu=length(unique(analdata$Mthd))
dimdir=16; fitcoefdir=matrix(0,nu/2,dimdir+1)
dimun=8; fitcoefun=matrix(0,nu/2,dimun+1)
j=0; k=0
for (i in 1:nu) {
 if (directed[i]==1) {
  j=j+1
  mod <- glm(ego_icap~size+nisolate+Ncomponent+diam+density_str+centralization+cor_inout+assort_inin+assort_inout+assort_outin+assort_outout+recip+trans+avg_clust+sumTrans+sum3Cycle,family=binomial(link="logit"),data=analdata,subset=(Mthd==i))
  #Augment with overall summary measure of fit (c-statistic)
  cstat <- OptimisedConc(mod)
  cprob <- as.numeric(cstat)[1]/100
  fitcoefdir[j,]=t(c(summary(mod)$coef[2:(dimdir+1),3],cprob))
 } else {
  k=k+1
  mod <- glm(ego_icap~size+nisolate+Ncomponent+diam+density_str+centralization+trans+avg_clust,family=binomial(link="logit"),data=analdata,subset=(Mthd==i))
  #Augment with overall summary measure of fit (c-statistic)
  cstat <- OptimisedConc(mod)
  cprob <- as.numeric(cstat)[1]/100
  fitcoefun[k,]=t(c(summary(mod)$coef[2:(dimun+1),3],cprob))
 }
}
fitcoefdir <- data.frame(fitcoefdir,row.names=NULL)
fitcoefun <- data.frame(fitcoefun,row.names=NULL)
names(fitcoefdir) <- c("size","nisolate","Ncomponent","diam","density_str","centralization","cor_inout","assort_inin","assort_inout","assort_outin","assort_outout","recip","trans","avg_clust","sumTrans","sum3Cycle","cprob")
names(fitcoefun) <- c("size","nisolate","Ncomponent","diam","density_str","centralization","trans","avg_clust","cprob")

#Now estimate models that control for the network predictors in the base model.
modcompdir=matrix(0,nu/2,7)
modcompun=matrix(0,nu/2,7)
indmissdir <- is.na(analdata$sizeb+analdata$nisolateb+analdata$Ncomponentb+analdata$diamb+analdata$density_strb+analdata$centralizationb+analdata$transb+analdata$avg_clustb+analdata$size+analdata$nisolate+analdata$Ncomponent+analdata$diam+analdata$density_str+analdata$centralization+analdata$cor_inout+analdata$assort_inin+analdata$assort_inout+analdata$assort_outin+analdata$assort_outout+analdata$recip+analdata$trans+analdata$avg_clust+analdata$sumTrans+analdata$sum3Cycle)
indmissun <- is.na(analdata$sizeb+analdata$nisolateb+analdata$Ncomponentb+analdata$diamb+analdata$density_strb+analdata$centralizationb+analdata$transb+analdata$avg_clustb+analdata$size+analdata$nisolate+analdata$Ncomponent+analdata$diam+analdata$density_str+analdata$centralization+analdata$trans+analdata$avg_clust)
j=0; k=0
for (i in 1:nu) {
 if (directed[i]==1) {
  j=j+1
  mod0 <- glm(ego_icap~sizeb+nisolateb+Ncomponentb+diamb+density_strb+centralizationb+transb+avg_clustb,family=binomial(link="logit"),data=analdata,subset=(Mthd==i & indmissdir==0))
  mod <- glm(ego_icap~sizeb+nisolateb+Ncomponentb+diamb+density_strb+centralizationb+transb+avg_clustb+size+nisolate+Ncomponent+diam+density_str+centralization+cor_inout+assort_inin+assort_inout+assort_outin+assort_outout+recip+trans+avg_clust+sumTrans+sum3Cycle,family=binomial(link="logit"),data=analdata,subset=(Mthd==i))
  Dev <- anova(mod0,mod) #Model comparison based on deviance; contains 4 statistics of interest
  #Augment with overall summary measure of fit (c-statistic)
  cstat0 <- OptimisedConc(mod0)
  cprob0 <- as.numeric(cstat0)[1]/100
  cstat <- OptimisedConc(mod)
  cprob <- as.numeric(cstat)[1]/100
  modcompdir[j,]=matrix(as.numeric(c(Dev[1,2],Dev[2,2:4],cprob0,cprob,cprob-cprob0)),nrow=1)
 } else {
  k=k+1
  mod0 <- glm(ego_icap~sizeb+nisolateb+Ncomponentb+diamb+density_strb+centralizationb+transb+avg_clustb,family=binomial(link="logit"),data=analdata,subset=(Mthd==i & indmissun==0))
  mod <- glm(ego_icap~sizeb+nisolateb+Ncomponentb+diamb+density_strb+centralizationb+transb+avg_clustb+size+nisolate+Ncomponent+diam+density_str+centralization+trans+avg_clust,family=binomial(link="logit"),data=analdata,subset=(Mthd==i))
  Dev <- anova(mod0,mod) #Model comparison based on deviance; contains 4 statistics of interest
  #Augment with overall summary measure of fit (c-statistic)
  cstat0 <- OptimisedConc(mod0)
  cprob0 <- as.numeric(cstat0)[1]/100
  cstat <- OptimisedConc(mod)
  cprob <- as.numeric(cstat)[1]/100
  modcompun[k,]=matrix(as.numeric(c(Dev[1,2],Dev[2,2:4],cprob0,cprob,cprob-cprob0)),nrow=1)
 }
}
modcompdir <- data.frame(modcompdir,row.names=NULL)
modcompun <- data.frame(modcompun,row.names=NULL)
names(modcompdir) <- c("DevBase","DevAll","DfDiff","DevDiff","cprobBase","cprobAll","cprobDiff")
names(modcompun) <- c("DevBase","DevAll","DfDiff","DevDiff","cprobBase","cprobAll","cprobDiff")

#Attached design information
DsgnMat=data.frame(Mthd=seq(1,80),Type=type,Interm=interm,Revisit=revisit,Directed=directed,Binary=binary)
ind=(directed==1)
modcompdir=data.frame(DsgnMat[ind,],modcompdir,row.names=NULL)
modcompun=data.frame(DsgnMat[!ind,],modcompun,row.names=NULL)
fitcoefdir=data.frame(DsgnMat[ind,],fitcoefdir,row.names=NULL)
fitcoefun=data.frame(DsgnMat[!ind,],fitcoefun,row.names=NULL)

#Output files (c = cross-sectional; l = lagged)
#All: Included number of components and diameter among undirected measures
write.csv(fitcoefdir,paste(outdir,"AdoptDirMeasures2011lnrAll.csv",sep=""),row.names=FALSE)
write.csv(fitcoefun,paste(outdir,"AdoptUndirMeasures2011lnrAll.csv",sep=""),row.names=FALSE)
write.csv(modcompdir,paste(outdir,"AdoptDirComp2011lnrAll.csv",sep=""),row.names=FALSE)
write.csv(modcompun,paste(outdir,"AdoptUndirComp2011lnrAll.csv",sep=""),row.names=FALSE)