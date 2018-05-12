#' Average Causal Effect with Acceptability Index
#'
#' Computes the Acceptability index as follows:
#' m(A) = (0.5)*(Al + Ar)
#' w(A) = (0.5)*(Ar - Al)
#' Ai = (m(B)-m(A))/(w(B)+w(A))
#' Sengupta et al. (2000). On comparing interval numbers. European Journal of Operational Research, 127(1)
#'
#' @param cpdag The CPDAG with the undefined causal structure 
#' @param cov.mat The covariance matrix of the set of variables V
#' @param x The variable x to determine ICE_x|y-hat and ICE_y|x-hat
#' @param y The variable y to determine ICE_x|y-hat and ICE_y|x-hat
#' 
#' @return A data frame with information of the computed ICE, the acceptability index and the estimated causal direction. 
#' 
#' 

ace_ai <- function (cpdag,cov.mat, x, y){
  
  ## Causal effect x->y
  #Ac.eff <- ida(x,y,cov(dat),getGraph(cpdag),method = "global")
  #Ac.eff <- idaFast(x,y,cov.mat,getGraph(cpdag))
  Ac.eff <- tryCatch(idaFast(x,y,cov.mat,getGraph(cpdag)),
     error = function(e) {0})
  #Ac.eff <- cbind(0,Ac.eff)
  
  #to.del <- which(Ac.eff == 0)
  #Ac.eff <- Ac.eff[-to.del]
  
  
  I <- range(Ac.eff)
  I <- signif(I, digits = 4)
  if (I[1]==I[2])
      I[1] <- 0
  Il <- I[1]
  Ir <- I[2]
  mA <- (0.5) * (Il + Ir)
  wA <- (0.5) * (Ir - Il)
  
  #Bc.eff <- ida(y,x,cov(dat), getGraph(cpdag), method = "global",type = gtype)
  #Bc.eff <- idaFast(y,x,cov.mat, getGraph(cpdag))
  Bc.eff <- tryCatch(idaFast(y,x,cov.mat,getGraph(cpdag)),
          error = function(e) {0})
  #Bc.eff <- cbind(0,Bc.eff)
  
  #to.del <- which(Bc.eff == 0)
  #Bc.eff <- Bc.eff[-to.del]
  
  I <- range(Bc.eff)
  I <- signif(I, digits = 4)

  if (I[1]==I[2])
    I[1] <- 0
  Il <- I[1]
  Ir <- I[2]
  mB <- (0.5) * (Il + Ir)
  wB <- (0.5) * (Ir - Il)
  
  stopifnot(!is.logical(mA))
  stopifnot(!is.logical(mB))
  if(is.infinite(mA) | is.infinite(mB))
    browser()
  if (mA > mB) {
    tmpm <- mA
    tmpw <- wA
    mA   <- mB
    wA   <- wB
    mB   <- tmpm
    wB   <- tmpw
    tmp  <- x
    x    <- y
    y    <- tmp
    tmpce <- Ac.eff
    Ac.eff <- Bc.eff
    Bc.eff <- tmpce
  }
  if((wB+wA)!=0){
  
    Ai <- (mB - mA) / (wB + wA)
    Ai <- signif(Ai, digits = 4)
    dir.est <- sprintf("%i,%i",y,x)
  }else{
    Ai <- 0
    dir.est <- sprintf("%i,%i",y,x)
  }

  #mat.interv <- matrix(c(x, y, mA, wA, y, x, mB, wB), ncol = 4, byrow = TRUE)
  #rownames(mat.interv)<-list(sprintf("%i->%i",x,y),sprintf("%i->%i",y,x))

  #dir.est (y,x) premise: y->x Ai corresponds to such premise.
  df.out <- data.frame(x=x,y=y,mA=mA,wA=wA,mB=mB,wB=wB,Ai=Ai,
                      dir.est=dir.est)
            #dir.real=NULL,ceff.real=NULL,p = NULL,pconn = NULL,conserv = NULL,nr.eUndir = NULL)
  
  return(df.out)
}