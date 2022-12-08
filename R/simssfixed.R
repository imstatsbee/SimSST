#
# ./R/Updated-SST-Simulation-Fixed SSDs.R
#
# Chel Hee Lee & Mohsen Soltanifar
# 2022-DEC-03
#

#' @rdname simssfixed
simssfixed0 <- function(pid, n, m, SSD.b, dist.go, theta.go,
                                      dist.stop, theta.stop){

  SRRT0 <- -999 # ???
  SSD1 <- SSD.b
  id1 <- as.vector(matrix(pid,nrow=1,ncol = n))

  if(!(dist.go %in% c("ExG", "SW"))) stop("This is incorrect Go distribution!")
  if(!(dist.stop %in% c("ExG", "SW"))) stop("This is incorrect Stop distribution!")

  if(dist.go == "ExG" & dist.stop == "ExG"){
    GORT1 <- round(gamlss.dist::rexGAUS(n, mu = theta.go[1], sigma = theta.go[2] , nu = theta.go[3]), digits = 1)
    SSRT1 <- round(gamlss.dist::rexGAUS(n, mu = theta.stop[1], sigma = theta.stop[2], nu = theta.stop[3]), digits = 1)
  }

  if(dist.go == "ExG" & dist.stop == "SW"){
    GORT1 <- round(gamlss.dist::rexGAUS(n, mu = theta.go[1], sigma = theta.go[2] , nu = theta.go[3]), digits = 1)
    SSRT1 <- round(gamlss.dist::rIG(n, mu = theta.stop[1], sigma = theta.stop[2])+ theta.stop[3], digits = 1)
  }

  if(dist.go == "SW" & dist.stop == "ExG"){
    GORT1 <- round(gamlss.dist::rIG(n, mu = theta.go[1], sigma = theta.go[2])+ theta.go[3], digits = 1)
    SSRT1 <- round(gamlss.dist::rexGAUS(n, mu = theta.stop[1], sigma = theta.stop[2], nu = theta.stop[3]), digits = 1)
  }

  if(dist.go == "SW" & dist.stop == "SW"){
    GORT1 <- round(gamlss.dist::rIG(n, mu = theta.go[1], sigma = theta.go[2])+ theta.go[3], digits = 1)
    SSRT1 <- round(gamlss.dist::rIG(n, mu = theta.stop[1], sigma = theta.stop[2])+ theta.stop[3], digits = 1)
  }


  SSRT1SSD1 <- SSRT1+SSD1
  SSD1 <- as.vector(matrix(SSD1, nrow=1, ncol=n))
  SRRT1 <- as.vector(matrix(SRRT0, nrow=1, ncol=n))

  TrialType1 <- as.vector(matrix('Go', nrow=1, ncol = n))
  Inhibition1 <- as.vector(matrix('-999', nrow=1, ncol = n))

  for (i in 1:m)
  {
    if (GORT1[i] > SSRT1SSD1[i]) TrialType1[i]='Stop(Successful)'
    else TrialType1[i]='Stop(Failed)'

    if (GORT1[i] > SSRT1SSD1[i]) Inhibition1[i]='1'
    else Inhibition1[i]='0'
  }

  for (i in (m+1):n)
  {
    if (SSRT1[i]!= 'NA') SSRT1[i]= -999
    if (SSRT1SSD1[i] != 'NA') SSRT1SSD1[i]= -999
    if (SSD1[i] != 'NA') SSD1[i]= -999
  }

  SRRT1 <- ifelse(TrialType1 %in% c('Stop(Failed)'), GORT1, SRRT0)
  GORT1 <- ifelse(TrialType1 %in% c('Stop(Successful)','Stop(Failed)'), -999, GORT1)

  TrialType2 <- dplyr::recode(TrialType1,
                              'Stop(Successful)' = "Stop",
                              'Stop(Failed)' = "Stop",
                              'Go' = "Go")

  MAT1 <- matrix(c(id1, TrialType2, Inhibition1, GORT1, SSRT1,SRRT1, SSD1), nrow=length(GORT1))
  MAT11 <- MAT1[sample(nrow(MAT1)),]

  colnames(MAT11) <- c('Participant.id', 'Trial','Inhibition', 'GORT', 'SSRT', 'SRRT','SSD')

  return(MAT11)
}

#' @rdname simssfixed
#' @title Simulatng SSRT data using fixed SSD methods
#' @description This function simulates b>=1 blocks of stop signal task trials for several participants using fixed SSD method.
#' @param pid a vector of size b of Participant.id
#' @param block a block name vector of size b blocks
#' @param n a vector of size b of total number of trials
#' @param m a vector of size b of total number of stops
#' @param SSD.b a vector of size b of stop signal delay
#' @param dist.go a vector of size b of distribution of go trials (ExG or SW)
#' @param dist.stop a vector of size b of distribution of stop.trials (ExG or SW)
#' @param theta.go c(mu.go, sigma.go, tau.go), a b*3 matrix
#' @param theta.stop c(mu.stop, sigma.stop, tau.stop), a b*3 matrix
#' @returns M11
#' Output: a giant matrix with "sum(n)" rows and (7+1) columns
#'
#' @references
#' Gordon D. Logan. On the Ability to Inhibit Thought and Action: A User's Guide to the Stop Signal Paradigm. In D. Dagenbach, & T.H. Carr (Eds.), Inhibitory Process in Attention, Memory and Language. San Diego: Academic Press, 1994.
#'
#' Dora Matzke, Jonathon Love, Thomas V. Wiecki, Scott D. Brown, and et al. Release the BEESTS: Bayesian Estimation of Ex-Gaussian Stop Signal Reaction Times Distributions. Frontiers in Psychology, 4: Article 918, 2013.
#'
#' Mohsen Soltanifar. Stop Signal Reaction Times: New Estimations with Longitudinal, Bayesian and Time  Series based Methods, PhD Dissertation, Biostatistics Division, Dalla Lana School of Public Health, University of Toronto, Toronto, Canada, 2020.
#'
#' @examples
#' mySSTdata1 <- simssfixed(
#'  pid = c("John.Smith","Jane.McDonald","Jane.McDonald"),
#'  n = c(50,100,150), m=c(10,20,30), SSD.b=c(200,220,240),
#'  dist.go=c("ExG","ExG","ExG"),
#'  theta.go=as.matrix(rbind(c(400,60,30),c(440,90,90),c(440,90,90))),
#'  dist.stop=c("ExG","ExG","ExG"),
#'  theta.stop=as.matrix(rbind(c(100,70,60),c(120,80,70),c(120,80,70))),
#'  block=c(1,1,2))
#' mySSTdata1
#'
#' @export

# SstSimulatedFixedSsd
simssfixed <- function(pid, block, n, m, SSD.b,
                                 dist.go, theta.go,
                                 dist.stop, theta.stop
                                 ) {

  b <- length(block)
  csn <- c(0, cumsum(n))
  M1 <- matrix(NA, nrow = sum(n), ncol = 8)

  for(i in 1:b){
    M1[c((csn[i]+1):csn[i+1]),1] <- block[i]
    M1[c((csn[i]+1):csn[i+1]),c(2:8)] <- simssfixed0(pid=pid[i], n=n[i], m=m[i], SSD.b=SSD.b[i], dist.go=dist.go[i], theta.go=theta.go[i,], dist.stop=dist.stop[i], theta.stop=theta.stop[i,])
  }

  M1[,c(1,2)] <- M1[,c(2,1)]
  M11 <- M1
  colnames(M11) <- c('Participant.id','Block','Trial','Inhibition','GORT','SSRT','SRRT','SSD')

  return(M11)
}

