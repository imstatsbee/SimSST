## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
rm(list=ls())
library(SimSST)
library(dplyr)
set.seed(1)

## -----------------------------------------------------------------------------
# pid <- c("FNLN1","FNLN1")
# block <- c(1,2)
# n <- c(10,10) 
# m <- c(4,4)   
# SSD.b <- c(220,240)
# dist.go <- c("ExG","ExG")
# theta.go <- as.matrix.data.frame(rbind(c(440,90,90),c(440,90,90))) 
# dist.stop <- c("ExG","ExG")
# theta.stop <- as.matrix.data.frame(rbind(c(120,80,70),c(120,80,70)))
mySSTdata1 <- 
  SstSimulatedFixedSsd(
    pid=c("FNLN1","FNLN1"),
    n=c(10,10), m=c(4,4), SSD.b=c(220,240),
    dist.go=c("ExG","ExG"),
    theta.go=as.matrix.data.frame(rbind(c(440,90,90),c(440,90,90))),
    dist.stop=c("ExG","ExG"),
    theta.stop=as.matrix.data.frame(rbind(c(120,80,70),c(120,80,70))), 
    block = c(1,2))
#View(mySSTdata1)
mySSTdata1 

## -----------------------------------------------------------------------------
# pid<-c("FNLN1","FNLN1")
# block<-c(1,2)
# n<-c(10,10) 
# m<-c(4,4)   
# SSD.b<-c(220,240)
# dist.go=c("ExG","ExG" )
# theta.go= as.matrix.data.frame(rbind(c(440,90,90),c(440,90,90))) 
# dist.stop=c("ExG","ExG" )
# theta.stop=as.matrix.data.frame(rbind(c(120,80,70),c(120,80,70)))

mySSTdata2 <- 
  SstSimulatedTrackingMethod(
    pid=c("FNLN1","FNLN1"),
    n=c(10,10), m=c(4,4), SSD.b=c(220,240), 
    dist.go=c("ExG","ExG" ),
    theta.go=as.matrix.data.frame(rbind(c(440,90,90),c(440,90,90))),
    dist.stop=c("ExG","ExG" ),
    theta.stop=as.matrix.data.frame(rbind(c(120,80,70),c(120,80,70))), 
    block=c(1,2))
# View(mySSTdata2)
mySSTdata2

## -----------------------------------------------------------------------------
Datatemp2<-mySSTdata2 
ss_presented<-recode(Datatemp2[,3], 'Stop' = "1", 'Go' = "0")
inhibited<-Datatemp2[,4]
ssd<-Datatemp2[,8]
rt<-Datatemp2[,5]
srrt<-Datatemp2[,7]
Data2<-cbind.data.frame(ss_presented,inhibited,ssd,rt,srrt)
for(i in 1:20){ if(Data2$inhibited[i]==0){Data2$rt[i]<-Data2$srrt[i]}}
myBEESTSdata2<-(Data2[,-5])[order(ss_presented),]
# View(myBEESTSdata2)
myBEESTSdata2

