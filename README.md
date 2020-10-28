# slice
Structure Learning via Intervals of Causal Effects (SLICE)

Resolve undefined causal edges froma an equivalence class represented bya a completed partially directed acyclic graphCPDAG
Params
 G.am An adjacency matrix of the initial CPDAG containig one or more undefined causal edges to be resolved. The edgemark-code in am refers to the row index.
 V.data A matrix with "m" observations (rows) of the "p" variables (columns)
Return
 A data frame with the learnt structure ("G.est"), the type of outcome ("G.outcome"), 
a data frame with information regarding the computed intervals of causal effects ("G.df_ice"),
the number of iteration elapsed ("iters").
Author: Samuel Montero-Hernandez, s.monterohdz@gmail.com


## Examples
source("slicegm.R") 
source("ace_ai-v2.R")
source("miscfunctions.R")

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
CPDAG.est.am <- as(CPDAG.est@@graph,"matrix")
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
