

setwd("~/Documents/WORK/IMAS/Sea_urchins/Bioeconomic_model/working model/Paper 1/Best/FitSTAN_Best")

#' @title STM - Generates the Size Transition Matrix for Inverse Logistic
#'
#' @description STM - With the input of the four parameters inside a vector, 
#'     and a vector of initial lengths or mid-points of size classes STM 
#'     generates a square transition matrix with the probabilities of 
#'     growing from each initial size into the same or larger sizes. 
#'     Negative growth is disallowed. All columns in the matrix sunm to one.
#' @param p a vector of four parameters in the following order
#'     MaxDL - the maximum growth increment of the inverse logistic,
#'     L50 - the initial length at which the inflexion of the growth 
#'     increment curve occurs, L95 - the initial length that defines the 
#'     95th percentile of the growth increments, SigMax - the maximum 
#'     standard deviaiton of the normal distribution used to describe the
#'     spread of each distribution of growth increments
#' @param mids - a vector of initial lengths which also define the width of 
#'     each size class thus, mids from 2 - 210 woul dimply 2mm size 
#'     classes 1 - 3 = 2, 3 - 5 = 4, etc
#' @return A square matrix with dimension =the length of the mids vector
#' @export 
#' 
#' @references Haddon, M., Mundy, C., and D. Tarbath (2008) Using an 
#'     inverse-logistic model to describe growth increments of blackip 
#'     abalone (Haliotis rubra) in Tasmania. Fisheries Bulletin 106: 58-71
#' @examples
#' \dontrun{
#' param <- c(25.0,120.0,170.0,4.0)
#' midpts <- seq(2,210,2)
#' G <- STM(param,midpts)
#' print(round(G[1:30,1:8],4))
#' }
#' 
#rm(list=ls())
#nmax <- 26
sizemax <- 150
sizemin <- 25
nintervals <- 50
#dev.off()
#MaxDL = 2.622
#L50 = 17.369
#L95 = 27.663
#stdevmax = 0.437
param <- c(2.622,17.369,27.633,1,.284,.451) ####################### Ling actual values
#param <- c(3,17.369,27.633,1,.284,.451) ####################### Ling actual values
#param <- c(2.622,17.369,30,1,.284,.451) ########################### new ones
#param <- c(3,17.369,27.633,1,.284,.451) # try MaxDL 3mm
#param <- c(3.5,17.369,27.633,1,.284,.451) # try MaxDL 3.5mm
#param <- c(4,17.369,27.633,1,.284,.451) # try MaxDL 4mm
#param <- c(2.88,19.1,30.423,1,.284,.451) # try MaxDL 4mm ### plus 10% to Ling
#param <- c(3.15,20.84,33.1596,1,.284,.451) # try MaxDL 4mm ### plus 20% to Ling
step <- (sizemax-sizemin)/nintervals
tdmids <- seq(sizemin+step/2,sizemax,step)
midpts <- tdmids/4.14
nmax <- length(midpts)

30*4.14
27.33*4.14
140/4.14
#pcheck <- param*1.2


# checking calculations
SigL = param[5]*midpts^param[6]
MeanL <- midpts + (param[1]/((1+exp(log(19.0)*(midpts-param[2])/(param[3]-param[2])))))



STM <- function(p,mids) { #    # p <- popparam[1:4]; mids <- midpts
  n <- length(mids)
  G <- matrix(0,nrow=n,ncol=n)
  cw <- mids[2]-mids[1]
  #SigL = p[5]*mids^p[6]
  growthincL  <- p[1]/((1+exp(log(19.0)*(mids-p[2])/(p[3]-p[2]))))
  SigL <- p[5]*growthincL^p[6]
  MeanL <- (mids + (p[1]/((1+exp(log(19.0)*(mids-p[2])/(p[3]-p[2]))))))
  for (j in 1:n) {
    for (i in 1:n) {
      Prob <- (1-pnorm(mids[i]+cw/2.0,MeanL[j],SigL[j],FALSE))
      if (i < j)  { G[i,j] <- 0.0 }
      if (i == j) { G[i,j] <- Prob }
      if (i > j)  { G[i,j] <- (Prob - (1-pnorm(mids[i-1]+cw/2.0,MeanL[j], SigL[j],FALSE))) }
    }
  }
  G[n,] <- G[n,]+ (1-colSums(G)) # plus group 
  rownames(G) <- mids
  colnames(G) <- mids
  class(G) <- "STM"
  return(G)
} # end of STM


G <- STM(param,midpts)
G[1:nmax,1:nmax]

str(G)
growthtransition <- as.data.frame(G[1:nmax,1:nmax])
sum(G[,16])

G[40:50,40:50]

#print(round(G[1:5,1:5],10))

Gparams <- param
Gparams[7:nmax] <- 0
Gparams <- rbind(Gparams,G[1:nmax,1:nmax])

#write.csv(Gparams,file=paste0(param[1],"MaxDL.csv"))
write.csv(G[1:nmax,1:nmax],file=paste0(param[1],"MaxDL_JL.csv"))

          