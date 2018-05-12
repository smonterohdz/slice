#' Structure Learning via Intervals of Causal Effects (SLICE)
#'
#' Resolve undefined causal edges froma an equivalence class represented bya a completed partially directed acyclic graphCPDAG
#' 
#' @param G.am An adjacency matrix of the initial CPDAG containig one or more undefined causal edges to be resolved. The edgemark-code in am refers to the row index.
#' @param V.data A matrix with \code{"m"} observations (rows) of the \code{"p"} variables (columns)
#'
#' @return A data frame with the learnt structure (\code{"G.est"}), the type of outcome (\code{"G.outcome"}), 
#' a data frame with information regarding the computed intervals of causal effects (\code{"G.df_ice"}),
#' the number of iteration elapsed (\code{"iters"}).
#' @author Samuel Montero-Hernandez, \email{samuel@@inaoep.mx}
#' @keywords graphical models
#'
#' @examples
#' source("slicegm.R") 
#' source("ace_ai-v2.R")
#' source("miscfunctions.R")
#' 
#' set.seed(101)
#' p <- 7
#' DAG.true <- randomDAG(p, prob = 0.5) 
#' CPDAG.true <- dag2cpdag(DAG.true) 
#' 
#' ## generate a set of observations
#' n.samp <- 10000
#' data.samp <- rmvDAG(n.samp, DAG.true)
#' 
#' ## estimate a CPDAG with PC algorithm
#' suffStat <- list(C = cor(data.samp), n = n.samp)
#' CPDAG.est <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p=p)
#' 
#' par(mfrow=c(1,3))
#' plot(DAG.true,main="TrueModel")
#' plot(CPDAG.true,main="CPDAGtrue")
#' plot(CPDAG.est,main="CPDAGestim")
#' 
#' 
#' ## Determine the causal edge direction from an estimated CPDAG
#' CPDAG.est.am <- as(CPDAG.est@@graph,"matrix")
#' slice.res <- slicegm(CPDAG.est.am,data.samp) 
#' par(mfrow=c(1,3))
#' plot(DAG.true,main="TrueModel")
#' plot(CPDAG.est,main="CPDAGestim")
#' plot(slice.res$G.est,main="SLICEoutcome")
#' shd.rates <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.est)
#' print(sprintf("From an estimated CPDAG: SHD=%.2f",shd.rates$shd.norm))
#' 
#' ## Determine the causal edge direction from a given true CPDAG
#' CPDAG.true.am <- as(CPDAG.true,"matrix")
#' slice.res <- slicegm(CPDAG.true.am,data.samp) 
#' par(mfrow=c(1,3))
#' plot(DAG.true,main="TrueModel")
#' plot(CPDAG.true,main="CPDAGtrue")
#' plot(slice.res$G.est,main="SLICEoutcome")
#' shd.rates <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.true)
#' print(sprintf("From a true CPDAG: SHD=%.2f",shd.rates$shd.norm))
#' 
#' @export
#' 




slicegm <- function(G.am, V.data) {
  amGM <- graphAM(adjMat = G.am,edgemode = "directed")
  G.cpdag <- as(amGM,"graphNEL")
  cpdagFlag <- FALSE
  beta <- (-1)
  iters <- 0
  
  df.null <- data.frame(
    x = NULL,
    y = NULL,
    mA = NULL,
    wA = NULL,
    mB = NULL,
    wB = NULL,
    Ai = NULL,
    dir.est = NULL,
    dir.real = NULL,
    ceff.real = NULL,
    p = NULL,
    pconn = NULL,
    conserv = NULL,
    nr.eUndir = NULL
  )
  
  # covariance matrix from data
  cov.mat = cov(V.data)
  
  # if (class(G.cpdag) == "graphNEL")
  #   am <- as(G.cpdag, "graphAM")@adjMat
  # if (class(G.cpdag) == "pcAlgo")
  #   am <- as(G.cpdag, "amat")
  
  ## list of undefined edges
  undir.res <- get.undiredges(G.am)
  e.list <- undir.res$e.list
  G.df_ice <- df.null
  nr.eUndir <- nrow(e.list)
  
  while (nr.eUndir > 0 & cpdagFlag == FALSE) {
    iters <- iters + 1
    df.iceL <- df.null
    for (i in 1:nr.eUndir) {
      ## Causal effect x->y
      x <- e.list[i, 1]
      y <- e.list[i, 2]
      
      tmpdf <- ace_ai(
        cpdag = G.cpdag,
        cov.mat = cov.mat,
        x = x,
        y = y
      )
      # 'df.aiceff' stores Ai, ceff.ral values (among others values) for every undefined edge of the current CPDAG
      df.iceL <- rbind(df.iceL, tmpdf)
      # 'df.aiceff' saves the relation between Ai and causal effects
      #df.aiceff <- rbind(df.aiceff,data.frame(Ai=tmpdf$Ai,c.eff=ceff.real,conserv=conserv))
    }# end undefined edges iteration
    G.df_ice <- rbind(G.df_ice,df.iceL)
    ## Orientation according maximum Ai or thresholding with beta
    ornt.res <- orientation(df.iceL, beta, G.cpdag)
    ## New undefined edges
    e.list <- ornt.res$e.list
    nr.eUndir <- nrow(e.list)
    cpdagFlag <- ornt.res$cpdagFlag
    G.cpdag <- ornt.res$cpdag.struct
  }# en while n.eundir >0 AND cpdagFlag==FALSE
  if (cpdagFlag == TRUE) {
    outcome <- "CPDAG"
  } else{
    outcome <- "DAG"
  }

  return(list(
    G.est = G.cpdag,
    G.outcome = outcome,
    G.df_ice = G.df_ice,
    iters = iters
  ))
  
}