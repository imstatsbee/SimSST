#
# ./R/Updated-SST-Simulation-TrackingMethod.R
#
# Chel Hee Lee & Mohsen Soltanifar
# 2022-DEC-03
#

#' @title One Block Model
#' @description Step(1): Do for one block and ExG/SW distribution
#' (for either Go/Stop)   OK
#' @details
#' @param pid Participant.id
#' @param n total number of trials (positive integer)
#' @param m total number of stops  (positive integer)
#' @param dist.go distribution of go trials    (ExG or SW)
#' @param dist.stop distribution of stop.trials  (ExG or SW)
#' @param theta.go c(mu.go,sigma.go,tau.go)
#' @param theta.stop c(mu.stop,sigma.stop,tau.stop)
#' @param SSD.b one block starting stop signal delay
#' @references
#' @author
#' @returns MAT22
#' @examples
#' ##
#' ## Example1
#' ##
#' myblockSSTdata1 <- SstSimulatedTrackingMethodBlock(
#'     pid="John.Smith", n=50, m=10, SSD.b=200, dist.go="ExG,
#'     theta.go=c(400,60,), dist.stop="ExG", theta.stop=c(100,70,60))
#' myblockSSTdata1
#'
#' ##
#' ## Example2
#' ##
#'
#' myblockSSTdata2 <- SstSimulatedTrackingMethodBlock(
#'      pid="Jane.Smith", n=50, m=10, SSD.b=200, dist.go="SW",
#'      theta.go=c(100,0.01,50), dist.stop="SW", theta.stop=c(75,0.01,100))
#' myblockSSTdata2
#'
#' ##
#' ## Example3
#' ## This example is expected to created an error
#' ##
#' myblockSSTdata3 <- SstSimulatedTrackingMethodBlock(
#' pid="Jane.Smith", n=50, m=10, SSD.b=200, dist.go="LN",
#' theta.go=c(100,0.01,50), dist.stop="SW", theta.stop=c(75,0.01,100))
#' myblockSSTdata3
#'
#'
#' @export

SstSimulatedTrackingMethodBlock <- function(pid, n, m, SSD.b,
                                            dist.go, theta.go,
                                            dist.stop, theta.stop)
{

  SRRT0 <- -999
  index <- as.vector(matrix(0,nrow=1,ncol = n))
  id <- as.vector(matrix(pid,nrow=1,ncol = n))

  if(!(dist.go %in% c("ExG", "SW"))) {
    warning("This is incorrect Go distribution!") & break
  }
  if(!(dist.stop %in% c("ExG", "SW"))) {
    warning("This is incorrect Stop distribution!") & break
  }

  if(dist.go=="ExG" & dist.stop=="ExG") {
    GORT  <- round(rexGAUS(n, mu = theta.go[1], sigma = theta.go[2] , nu = theta.go[3]), digits = 1)
    SSRT  <- round(rexGAUS(n, mu = theta.stop[1], sigma = theta.stop[2], nu = theta.stop[3]), digits = 1)
  }

  if(dist.go=="ExG" & dist.stop=="SW") {
    GORT  <- round(rexGAUS(n, mu = theta.go[1], sigma = theta.go[2] , nu = theta.go[3]), digits = 1)
    SSRT  <- round(rIG(n, mu = theta.stop[1], sigma = theta.stop[2])+ theta.stop[3], digits = 1)
  }

  if(dist.go=="SW" & dist.stop=="ExG") {
    GORT  <- round(rIG(n, mu = theta.go[1], sigma = theta.go[2])+ theta.go[3], digits = 1)
    SSRT  <- round(rexGAUS(n, mu = theta.stop[1], sigma = theta.stop[2], nu = theta.stop[3]), digits = 1)
  }

  if(dist.go=="SW" & dist.stop=="SW") {
    GORT  <- round(rIG(n, mu = theta.go[1], sigma = theta.go[2])+ theta.go[3], digits = 1)
    SSRT  <- round(rIG(n, mu = theta.stop[1], sigma = theta.stop[2])+ theta.stop[3], digits = 1)
  }

  SSD <- as.vector(matrix(0 ,nrow=1,ncol = n))
  SSD[1] <- SSD.b
  SRRT <- as.vector(matrix(SRRT0,nrow=1,ncol=n))

  TrialType <- as.vector(matrix('Go',nrow=1,ncol = n))
  Inhibition <- as.vector(matrix('-999',nrow=1,ncol = n))

  ##
  ## First: Check the fist stop trial
  ##

  if (GORT[1] > SSRT[1]+SSD[1]) {
    TrialType[1]<- 'Stop(Successful)'
  } else  {
    TrialType[1]<- 'Stop(Failed)'
  }
  if (GORT[1] > SSRT[1]+SSD[1]) {
    Inhibition[1]<- '1'
  } else  {
    Inhibition[1]<- '0'
  }
  if (GORT[1] > SSRT[1]+SSD[1]) {
    index[1]<-  +50
  } else  {
    index[1]<- -50
  }
  if (GORT[1] > SSRT[1]+SSD[1]) {
    SSD[2]<-  SSD[1]+50
  } else  {
    SSD[2]<- SSD[1]-50
  }

  ##
  ## Second: Use information above to create a martingale of SSDs !
  ##

  for (i in 1:(m-1)) {
    if (index[i]== 50)  {
      SSD[i+1]<- SSD[i]+50
    }
    if (index[i]== 50)  {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) {
        TrialType[i+1]<- 'Stop(Successful)'
      } else  {
        TrialType[i+1]<- 'Stop(Failed)'
      }
    }
    if (index[i]== 50)  {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) {
        index[i+1]<- +50
      } else {
        index[i+1]<- -50
      }
    }

    if (index[i]== -50) {
      SSD[i+1]<- SSD[i]-50
    }
    if (index[i]== -50) {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) {
        TrialType[i+1]<- 'Stop(Successful)'
      } else  {
        TrialType[i+1]<- 'Stop(Failed)'
      }
    }
    if (index[i]== -50) {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) {
        index[i+1]<- +50
      } else  {
        index[i+1]<- -50
      }
    }

    if (index[i]== 50)  {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) {
        Inhibition[i+1]<- '1'
      } else  {
        Inhibition[i+1]<- '0'
      }
    }
    if (index[i]== -50) {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) {
        Inhibition[i+1]<- '1'
      } else  {
        Inhibition[i+1]<- '0'
      }
    }
  }

  ##
  ## Third: Update Final Values of SSD given values of Index !
  ##

  SSRTSSD <- as.vector(matrix(0,nrow=1,ncol =n))

  for (i in 1:n) {
    SSRTSSD[i] <- SSRT[i] + SSD[i]
  }

  ##
  ## Fourth: We clear unwanted data for Go trials !
  ##

  for (i in (m+1):n) {
    if (SSRT[i]!= 'NA') {
      SSRT[i]= -999
    }
    if (SSRTSSD[i]!= 'NA') {
      SSRTSSD[i]= -999
    }
    if (SSD[i] != 'NA')  {
      SSD[i]= -999
    }
    if (index[i] != 'NA') {
      index[i]= -999
    }
  }

  SRRT <- ifelse(TrialType %in% c('Stop(Failed)'), GORT, SRRT0)
  GORT <- ifelse(TrialType %in% c('Stop(Successful)','Stop(Failed)'), -999, GORT)

  TrialTypee <- recode(TrialType, 'Stop(Successful)' = "Stop", 'Stop(Failed)' = "Stop", 'Go' = "Go")
  MAT2 <- matrix(c(id , TrialTypee, Inhibition , GORT , SSRT ,SRRT , SSD ), nrow=length(GORT ))

  MAT22 <- matrix(NA,ncol = 7, nrow=length(GORT ))
  for (i in 1:m){
    MAT22[2*i-1,] <- MAT2[i,];
    MAT22[2*i,] <- MAT2[m+i,]
  }
  for(i in (2*m+1):n) {
    MAT22[i,] <- MAT2[i,]
  }

  colnames(MAT22) <- c( 'Participant.id', 'Trial','Inhibition', 'GORT', 'SSRT', 'SRRT','SSD');
  return(MAT22)
}


#' @title B Blocks Model
#' @description Step(2): Do for b blocks and ExG/SW distribution
#' (for either Go/Stop)   OK
#' @details
#' @param block a block name vector of size b blocks
#' @param pid a vector of size b of Participant.id
#' @param n a vector of size b of total number of trials
#' @param m a vector of size b of total number of stops
#' @param dist.go  a vector of size b of   distribution of go trials    (ExG or SW)
#' @param dist.stop a vector of size b of   distribution of stop.trials  (ExG or SW)
#' @param theta.go c(mu.go,sigma.go,tau.go)   a b*3 matrix
#' @param theta.stop c(mu.stop,sigma.stop,tau.stop)   a b*3 matrix
#' @param SSD.b  a vector of size b of starting stop signal delay
#' @author
#' @references
#' @returns Output a giant matrix with "sum(n)" rows and (7+1) columns
#' @examples
#' #Example1
#' block<-c(1,1,2)
#' # FIXME:: this block parameter is not defined in this function
#' mySSTdata1<-SstSimulatedTrackingMethod(
#'     pid=c("John.Smith","Jane.McDonald","Jane.McDonald"),
#'     n=c(50,100,150), m=c(10,20,30),
#'     SSD.b=c(200,220,240), dist.go=c("ExG","ExG","ExG"),
#'     theta.go=as.matrix.data.frame(rbind(c(400,60,30),c(440,90,90),c(440,90,90))),
#'     dist.stop=c("ExG","ExG","ExG"),
#'     theta.stop=as.matrix.data.frame(rbind(c(100,70,60),c(120,80,70),c(120,80,70))))
#' mySSTdata1
#'
#' #Example2
#' block<-c(1,1,2)
#' mySSTdata2 <- SstSimulatedTrackingMethod(
#'   pid=c("John.Smith","Jane.McDonald","Jane.McDonald"),
#'   n=c(50,100,150), m=c(10,20,30), SSD.b = c(200,220,240),
#'   dist.go=c("ExG","ExG","ExG"),
#'   theta.go=as.matrix.data.frame(rbind(c(400,60,30),c(440,90,90),c(440,90,90))),
#'   dist.stop=c("SW","SW","SW"),
#'   theta.stop=as.matrix.data.frame(rbind(c(75,0.01,100),c(75,0.01,100),c(75,0.01,100))))
#' mySSTdata2
#'
#' # Example3  <--- produce error !
#'
#' mySSTdata3 <- SstSimulatedTrackingMethod(pid="John.Smith", n=c(1), m=c(50),
#'    SSD.b=c(200), dist.go=c("LN"),
#'    theta.go=as.matrix.data.frame(rbind(c(400,60,30))),
#'    dist.stop=c("SW"),
#'    theta.stop=as.matrix.data.frame(rbind(c(75,0.01,100))))
#' mySSTdata3
#'
#'
#' @export
#'

SstSimulatedTrackingMethod <- function(pid, n, m, SSD.b, dist.go, theta.go,
                                       dist.stop, theta.stop) {

  b<-length(block)
  csn<-c(0,cumsum(n))
  M2 = matrix(NA, nrow = sum(n), ncol = 8)

  for(i in 1:b){
    M2[c((csn[i]+1):csn[i+1]),1] <- block[i]
    M2[c((csn[i]+1):csn[i+1]),c(2:8)] <- SstSimulatedTrackingMethodBlock(pid[i],n[i],m[i],SSD.b[i],dist.go[i],theta.go[i,],dist.stop[i],theta.stop[i,])
  }

  M2[,c(1,2)] <- M2[,c(2,1)]
  M22 <- M2
  colnames(M22) <- c('Participant.id','Block','Trial','Inhibition','GORT','SSRT','SRRT','SSD')

  return(M22)
}

