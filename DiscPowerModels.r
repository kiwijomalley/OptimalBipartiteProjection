#Reads in data sets containing network-level statistics for each of 80 ways of building physician
# network for each PHN

library(lme4)
setwd("/projects/active/14593/programs/jomalley/SmartNetwork")
datdir="/projects/active/14593/idata/core_c/core_c_jomalley/PHN-subnetwork/2010/output-"
outdir="/projects/active/14593/idata/core_c/core_c_jomalley/PHN-subnetwork/2010/SmartNetwork/"

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

#Set missing values?
ind=(netdata==-99)
netdata[ind]=NA
indrow=(netdata$size<=20) #Due to Chuankai's code not computing measures when network has size<20 nodes
netdata[indrow,3:(ncol(netdata)-1)]=NA #Last column is netID; don't want it to be NA

#Add design variables based on Chuankai's scheme (easy to do here as just repeat for each of the Nnetworks)
# Will need to define these factors directly from Mthd to build into physician-level data later
netdata$Mthd=rep(seq(1,80),times=Nnetwork) #Because start count at 1 don't need to use sMthd==(i-1) later!
netdata$Type=rep(rep(seq(1,5),times=80/5),times=Nnetwork)
netdata$Interm=rep(rep(rep(c(0,1),each=5),times=80/10),times=Nnetwork) #Interm = 0 means intermediates not allowed
netdata$Revisit=rep(rep(rep(c(0,1),each=10),times=80/20),times=Nnetwork) #Revisit = 1 means revisit required
netdata$Directed=rep(rep(rep(c(1,0),each=20),times=80/40),times=Nnetwork)
netdata$Binary=rep(rep(rep(c(0,1),each=40),times=80/80),times=Nnetwork)

#Statistical transformations that help with comparability of the same measure across the 80 designs
totnodes <- netdata$size+netdata$nisolate #Because nisolate > size for 12.1% of observations, assume size refers to largest connected component
netdata$tsize <- totnodes #Don't overwrite size as totsize should be invariant across design methods but size of largest connected component won't necessarily be
netdata$nisolate <- netdata$nisolate/totnodes; names(netdata)[2] <- 'pisolate'
tottriads <- totnodes*(totnodes-1)*(totnodes-2)/6
alltriads <- apply(netdata[,13:27],1,'sum') #Sum's to 100000 as samples were taken as opposed to evaluating over the entire network
netdata$sumTrans <- netdata$sumTrans/100000 #Not perfect; get some values > 1
netdata$sum3Cycle <- netdata$sum3Cycle/100000 #Not perfect; get some values > 1
netdata$centralization <- netdata$centralization/((netdata$density_str*(netdata$size-1))^2) #Relative to degree

relwtsbase <- 66 #Base = binary, undirected, type 1, intermediate allowed and revisit not required
relwts <- tapply(netdata$density_str,netdata$Mthd,'mean',na.rm=TRUE)
relwtsvec <- rep(relwts/relwts[relwtsbase],times=nrow(netdata)/80) #80 = Number of distinct designs
netdata$density_str=netdata$density_str/relwtsvec #Makes density/strength have same mean as base design 
tapply(netdata$density_str,netdata$Mthd,'mean',na.rm=TRUE) #Check
tapply(netdata$centralization,netdata$Mthd,'mean',na.rm=TRUE) #Check

#Standardize network measures to make them comparable across measures (for pooled across measures analysis) by having SD of 1 across entire data set
stdev <- sqrt(apply(netdata[,1:31],2,'var',na.rm=TRUE))
ones <- rep(1,nrow(netdata))
scalemat <- matrix(ones,ncol=1) %*% t(matrix(stdev,ncol=1))
netdata[,1:31] <- netdata[,1:31]/scalemat

#Trim data to just hospitals for which they have non-missing data across all designs (get from Directed=1, Revisit=1, and Interm=0)
missdsgn <- (netdata$Directed==1 & netdata$Revisit==1 & netdata$Interm==0 & netdata$Binary==0 & netdata$Type==1) #The design most prone to missing values
indmissdir <- is.na(netdata$size+netdata$pisolate+netdata$density_str+netdata$centralization+netdata$cor_inout+netdata$assort_inin+netdata$assort_inout+netdata$assort_outin+netdata$assort_outout+netdata$recip+netdata$trans+netdata$avg_clust+netdata$sumTrans+netdata$sum3Cycle)
goodhosp <- unique(netdata$netid[(missdsgn & !indmissdir)])
ind <- (netdata$netid %in% goodhosp)
netdata <- netdata[ind,] #Ensures that use same set of hospitals for all analyses


#Make repeated measures dataset with the network measure being a predictor variable.
nmeasure=16
sPhnID=rep(netdata$netid,times=nmeasure)
sMthd=rep(netdata$Mthd,times=nmeasure)
sType=rep(netdata$Type,times=nmeasure)
sInterm=rep(netdata$Interm,times=nmeasure)
sRevisit=rep(netdata$Revisit,times=nmeasure)
sDirected=rep(netdata$Directed,times=nmeasure)
sBinary=rep(netdata$Binary,times=nmeasure)
pweight=netdata$tsize*nrow(netdata)/sum(netdata$tsize) #Keep sample-size on par with unweighted
#pweight[indrow]=NA
precweight=rep(pweight,times=nmeasure) #Use total size as most measures involve all nodes but make NA if size<=20
MeasureID=rep(seq(1,nmeasure),each=nrow(netdata))
MeasureVal=c(netdata$size,netdata$pisolate,netdata$density_str,netdata$centralization,netdata$cor_inout,netdata$assort_inin,netdata$assort_inout,netdata$assort_outin,netdata$assort_outout,netdata$recip,netdata$trans,netdata$avg_clust,netdata$Ncomponent,netdata$diam,netdata$sumTrans,netdata$sum3Cycle)
snetdata=data.frame(cbind(sPhnID,sMthd,sType,sInterm,sRevisit,sDirected,sBinary,precweight,MeasureID,MeasureVal))
ind=(snetdata$MeasureID==1 | snetdata$MeasureID==13 | snetdata$MeasureID==14)
#snetdata$precweight[ind]=1/snetdata$precweight[ind] #Invert weights for size of largest connected component
#ind=(snetdata$MeasureID==13 | snetdata$MeasureID==14)
#snetdata$precweight[ind]=1
snetdata$precweight=1

# Compare across multiple networks which designs or design factors yield networks that vary the most between
#  networks (of phns) being compared !!!

#Compute network with the greatest ICC for PHN with measures repeated
undir=c(1:4,11:14)
dir=c(5:10,15:16)
snetdatadir=snetdata[(snetdata$MeasureID %in% dir),]
nu=length(unique(snetdata$sMthd))
iccnetdata=matrix(0,nu,4)
iccnetdataund=matrix(0,nu,4)
iccnetdatadir=matrix(0,nu,4)
iccnetmeadata=matrix(0,nu*nmeasure,4)
for (i in 1:nu) {
  tmp=as.data.frame(VarCorr(lmer(MeasureVal~factor(MeasureID)+(1|sPhnID),weights=precweight,subset=(sMthd==i),data=snetdata,control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))),comp="Variance"))
  iccnetdata[i,]=c(i,tmp$vcov[1],tmp$vcov[2],tmp$vcov[1]/(tmp$vcov[1]+tmp$vcov[2]))
  snetdatasub1=snetdata[(snetdata$sMthd==i & snetdata$MeasureID %in% undir),]
  if (nrow(snetdatasub1)>sum(is.na(snetdatasub1$MeasureVal))) {
    tmp=as.data.frame(VarCorr(lmer(MeasureVal~factor(MeasureID)+(1|sPhnID),weights=precweight,data=snetdatasub1),comp="Variance"))
    iccnetdataund[i,]=c(i,tmp$vcov[1],tmp$vcov[2],tmp$vcov[1]/(tmp$vcov[1]+tmp$vcov[2]))
  } else {
    iccnetdataund[i,]=c(i,NA,NA,NA)
  }
  snetdatasub2=snetdata[(snetdata$sMthd==i & snetdata$MeasureID %in% dir),]
  if (nrow(snetdatasub2)>sum(is.na(snetdatasub2$MeasureVal))) {
    tmp=as.data.frame(VarCorr(lmer(MeasureVal~factor(MeasureID)+(1|sPhnID),weights=precweight,data=snetdatasub2),comp="Variance"))
    iccnetdatadir[i,]=c(i,tmp$vcov[1],tmp$vcov[2],tmp$vcov[1]/(tmp$vcov[1]+tmp$vcov[2]))
  } else {
    iccnetdatadir[i,]=c(i,NA,NA,NA)
  }
  for (j in 1:nmeasure) {
    ave=mean(snetdata$MeasureVal[sMthd==i & snetdata$MeasureID==j],na.rm=TRUE)
    tmp=var(snetdata$MeasureVal[sMthd==i & snetdata$MeasureID==j],na.rm=TRUE)
    pos=(i-1)*nmeasure+j
    iccnetmeadata[pos,]=c(i,j,ave,tmp)
  }
}
#Could alternatively control for sInterm, sRevisit, sDirected, sBinary
iccnetdata=data.frame(iccnetdata)
names(iccnetdata)=c("Mthd","bvar","wvar","icc")
write.csv(iccnetdata,paste(outdir,"ICCnetPHN2010nr.csv",sep=""),row.names=FALSE)
iccnetdataund=data.frame(iccnetdataund)
names(iccnetdataund)=c("Mthd","bvar","wvar","icc")
write.csv(iccnetdataund,paste(outdir,"ICCnetPHNund2010nr.csv",sep=""),row.names=FALSE)
iccnetdatadir=data.frame(iccnetdatadir)
names(iccnetdatadir)=c("Mthd","bvar","wvar","icc")
write.csv(iccnetdatadir,paste(outdir,"ICCnetPHNdir2010nr.csv",sep=""),row.names=FALSE)
iccnetmeadata=data.frame(iccnetmeadata)
names(iccnetmeadata)=c("Mthd","Measure","ave","wvar")
write.csv(iccnetmeadata,paste(outdir,"ICCnetmeaPHN2010nr.csv",sep=""),row.names=FALSE)