# The input adjacency matrix must has the format (the edgemark-code refers to the row index)
# 0, 1 for tail and arrowhead 
# 
# amat[a,b]=0 and amat[b,a]=1 implies a --> b
# amat[a,b]=1 and amat[b,a]=0 implies a <-- b
# amat[a,b]=0 and amat[b,a]=0 implies a     b
# amat[a,b]=1 and amat[b,a]=1 implies a <-> b


rm(list = ls())
library(pcalg)
library(Rgraphviz)

source("R/slicegm.R") 
source("R/ace_ai.R")
source("R/miscfunctions.R")

 set.seed(101)
 p <- 7
 DAG.true <- randomDAG(p, prob = 0.5) 
 CPDAG.true <- dag2cpdag(DAG.true) 
 
 ## generate a set of observations
 n.samp <- 10000
 data.samp <- rmvDAG(n.samp, DAG.true)
 
 ## estimate a CPDAG with PC algorithm
 suffStat <- list(C = cor(data.samp), n = n.samp)
 CPDAG.est <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p=p)
 
 par(mfrow=c(1,3))
 plot(DAG.true,main="TrueModel")
 plot(CPDAG.true,main="CPDAGtrue")
 plot(CPDAG.est,main="CPDAGestim")
 
 
 ## Determine the causal edge direction from an estimated CPDAG
 CPDAG.est.am <- as(CPDAG.est@graph,"matrix")
 slice.res <- slicegm(CPDAG.est.am,data.samp) 
 par(mfrow=c(1,3))
 plot(DAG.true,main="TrueModel")
 plot(CPDAG.est,main="CPDAGestim")
 plot(slice.res$G.est,main="SLICEoutcome")
 shd.rates <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.est)
 print(sprintf("From an estimated CPDAG: SHD=%.2f",shd.rates$shd.norm))
 
## Determine the causal edge direction from a given true CPDAG
 CPDAG.true.am <- as(CPDAG.true,"matrix")
 slice.res <- slicegm(CPDAG.true.am,data.samp) 
 par(mfrow=c(1,3))
 plot(DAG.true,main="TrueModel")
 plot(CPDAG.true,main="CPDAGtrue")
 plot(slice.res$G.est,main="SLICEoutcome")
 shd.rates <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.true)
 print(sprintf("From a true CPDAG: SHD=%.2f",shd.rates$shd.norm))
