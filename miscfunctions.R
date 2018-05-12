#============================================================================== 
## Generates a synthetic model with 'p' variables
synt.mdl <- function(p, pconn){
  zerodegFlag <- TRUE
  
  while (zerodegFlag ) {
    ## true DAG
    myDAG <- randomDAG(p, prob = pconn) 
    am <- (as(myDAG,"matrix")!=0)*1L
    
    if (any(colSums(am + t(am)) == 0)) {
      zerodegFlag = TRUE
      next
    } else{
      zerodegFlag = FALSE
    }
    
  }
  return(mdl=myDAG)
}

#============================================================================== 
## Generates a vector with non-gaussian error terms as in LiNGAM paper
# Shimizu et al. JMLR vol. 7, pp. 2003–2030, 2006
genErrDist <- function(p,nsamp,errminmax){
  # First, generate errstd
  #errstd = rand(dims,1)*(errminmax(2)-errminmax(1)) + errminmax(1);
  errstd <- array(runif(p) * (errminmax[2]-errminmax[1]) + errminmax[1],dim = c(p,1))
  
  # Nonlinearity exponent, selected to lie in [0.5, 0.8] or [1.2, 2.0].
  # (<1 gives subgaussian, >1 gives supergaussian)
  #q = rand(dims,1)*1.1+0.5;    
  #ind = find(q>0.8);           
  #q(ind) = q(ind)+0.4;     
  
  q <- matrix(runif(p) * 1.1 + 0.5,nrow = p)
  ind <- which(q>0.8,arr.ind = TRUE)
  q[ind] <- q[ind] + 0.4
  
  
  # Number of data vectors
  #samples = 10000;  
  nsamp
  
  # This generates the disturbance variables, which are mutually 
  # independent, and non-gaussian
  #S = randn(dims,nsamp);
  #S = sign(S).*(abs(S).^(q*ones(1,nsamp)));
  S <- matrix(rnorm(p * nsamp), nrow = p)
  S <- sign(S) * ( abs(S)**(q%*%array(1,nsamp))   )   
  
  # This normalizes the disturbance variables to have the 
  # appropriate scales
  #S = S./((sqrt( mean((S').^2)' )./disturbancestd)*ones(1,samples));
  S <- S / ( (sqrt(array(rowMeans(S**2),dim = c(p,1))) / errstd) %*% array(1,nsamp)   )
  return(t(S))
}


#============================================================================== 
## Generates a vector with non-gaussian error terms (error distrubances)
disturb <- function(type, p, nsamp, p1, p2){
  errMat <- switch(type,
    normal = matrix(rnorm(nsamp * p, mean = p1,sd = p2), nrow = nsamp),
    cauchy = matrix(rcauchy(nsamp * p, location = p1, scale = p2), nrow = nsamp),
    subsupergaussian = genErrDist(p = p,nsamp = nsamp,errminmax =c(0.5, 1.5))
  )
  return(errMat)
}

#==============================================================================
## Determine the type of a graph 'cpdag', 'pdag', 'dag'
which.graph <- function(am){
  require(matrixcalc)
  stopifnot(is.matrix(am),is.square.matrix(am))
  gtype<-NULL
  if (isValidGraph(am, "dag") == TRUE) {
    gtype <- "dag"
  } else{
    if (isValidGraph(am, "pdag") == TRUE){
      gtype <- "pdag"
    }else{
      if(isValidGraph(am, "cpdag") == TRUE){
        gtype <- "cpdag"
      }
    }
  }
  return(gtype)
}


#==============================================================================
## Determining undefined edges
get.undiredges <- function(am){
  #wm <- wgtMatrix(getGraph(pc.fit))
  p <- nrow(am)
  wmU <- am + t(am)
  e.list <- which(wmU == 2 & upper.tri(wmU), arr.ind = TRUE)
  nr.eUndir <- nrow(e.list)
  if (nr.eUndir == 0) {
    zeroundirFlag = TRUE
    #print("One memeber in the equivalence class!")
  } else{
    zeroundirFlag = FALSE
    #cat(nrow(e.list),"undefined edges to resolve in",p,"variables \n")
  }
  undir.res <- list(zeroundirFlag = zeroundirFlag, e.list=e.list)
  return(undir.res)
}


#==============================================================================
# Compares the true graph with an estimated graph in terms of 
# True Positive Rate (TPR), False Positive Rate (FPR) and True Discovery Rate (TDR).
mycomparegraphs <- function(amEst,amTrue){
  if(any(amEst==3)){
    amEst[which((amEst + t(amEst))<5)] <- 0
    amEst <- (amEst!=0)*1L
  }else{
    amEst[which((amEst + t(amEst))==2)] <- 0
  }
  
  TP <- sum((amEst & amTrue)*1L)
  FP <- sum( ((amEst - amTrue)==1)*1L )
  FN <- sum( ((amEst - amTrue)==-1)*1L )
  FE <- sum(amTrue)

  amEst  <- amEst  + t(amEst)
  amTrue <- amTrue + t(amTrue)
  diag(amEst)<-1
  diag(amTrue)<-1
  amEst  <- (!(amEst))*1L
  amTrue <- (!(amTrue))*1L
  TN <- sum((amEst & amTrue)*1L)/2
  TPR <- TP / (TP+FN)
  FPR <- FP / (FP+TN)
  TDR <- TP / FE
  if(is.nan(FPR)){
    FPR<- 0.0
    #browser()
  }
  if(is.nan(TDR)){
    TDR <- 0.0
  }
  
  return(data.frame(FE,TP,TN,FP,FN,TPR,FPR,TDR))
}


#==============================================================================
#' Function mycomparegraphs2 obtain C(ompleteness) and S(oundness)
#' #DefinedEdges = DirectedEdges_amEst- DirectedEdges_amCPDAG OR dirgidos de los #UndefinedEdges_amCPDAG 
#' #TotalEdges = DirectedEdges_amTrue - DirectedEdges_amCPDAG OR #UndefinedEdges_amCPDAG
#' #CorrectEdges = de los #DefinedEdges aquellos correctos (coinciden con el analogo en el modelo original)
#' C=#DefinedEdges/#TotalEdges and S=#CorrectEdges/#DefinedEdges
# edges codes
# class pcAlgo, as returned from skeleton() or pc()
# 0, 1 for (tail|noEdge) and arrowhead e.g. 
# Note that the edgemark-code refers to the row index (as opposed adjacency matrices of type mag or pag)
# amat[a,b]=0 and amat[b,a]=1 implies a --> b
# amat[a,b]=1 and amat[b,a]=0 implies a <-- b
# amat[a,b]=1 and amat[b,a]=1 implies a <-> b
# amat[a,b]=0 and amat[b,a]=1 implies a     b
# amat[b,a]=1 if a-->b

# class "fciAlgo" returned from fci(), rfci(), fciPlus(), and dag2pag(), 
# class "graphNEL" returned from addBgKnowledge(), addEdge(), removeEdge()
# class"LINGAM" returned from lingam()
#0,1,2,3 for no edge, circle, arrowhead, tail; e.g., 
# Note that the edgemark-code refers to the column index (as opposed adjacency matrices of type dag or cpdag).
# amat[a,b]=2 and amat[b,a]=3 implies a --> b
# amat[a,b]=3 and amat[b,a]=2 implies a <-- b
# amat[a,b]=2 and amat[b,a]=2 implies a <-> b
# amat[a,b]=1 and amat[b,a]=3 implies a --o b
# amat[a,b]=0 and amat[b,a]=0 implies a     b
# amat[a,b]=2 if a?->b
mycomparegraphs2 <- function(mdl.true,mdl.est,mdl.cpdag){
  #make amTrue <- wgtMatrix(myDAG,transpose = TRUE)
  #Determine the type of the graph
  amTrue  <- as(mdl.true,"graphAM")@adjMat
  C<-0
  S<-0
  if(is.null(mdl.cpdag)){
    amEst <- as(getGraph(mdl.est),"graphAM")@adjMat
    tot.edges <- sum(amTrue)
    def.edges <- sum(amEst)
    hit.edges <- sum(((amEst+ amTrue)==2)*1L)
  }else{
    amEst   <- as(getGraph(mdl.est),"graphAM")@adjMat
    amCPDAG <- as(getGraph(mdl.cpdag),"graphAM")@adjMat
    marks <- unique(as.vector(amEst))
    if(identical(marks,c(0,1))){
      typeEst<-"cpdag"
    }else{
      if(identical(marks,c(0,1,2,3))){
        typeEst <- "pag"
      }else{
        if(identical(marks,c(0))){
          typeEst <- "skel"
        }
      }
    }
    #TotalEdges from true model (assumes a DAG)
    uEdges <- get.undiredges(amCPDAG)
    e.list <- uEdges$e.list
    tot.edges <- nrow(e.list)
    #DefinedEdges
    eb.list <- rbind(e.list,e.list[,c(2,1)]) #edgeBidirected.list
    mm<-matrix(FALSE,nrow = nrow(amCPDAG),ncol = ncol(amCPDAG) ) #mask for edges of interest (bidirected)
    mm[eb.list] <- TRUE
    #
    if(typeEst=="skel"){
      return(list(C=0,S=0))
    }
    #DefinedEdges
    if(typeEst=="cpdag"){
      # only preserve the slots of bidirected edges in both amCPDAG and amEst
      amCPDAG[!(mm)] <- 0 
      amEst[!(mm)] <- 0 
      e.dir <- which((amEst+ t(amEst))==1,arr.ind = TRUE)
      def.edges <-nrow(e.dir)/2
      marks.est <- amEst[e.dir]
      marks.true <- amTrue[e.dir]
      hit.edges <- sum((marks.est==marks.true)*1L)/2
    }
  }
  if(tot.edges!=0){
    C <- def.edges/tot.edges
  }else{
    C <- 0
  }
  if(def.edges!=0){
    S <- hit.edges/def.edges
  }else{
    S <- 0
  }
  #if(beta==1 & nsamp==10000)
  #  browser()
  return(list(C=C,S=S))
}


#==============================================================================
# graph comparison based on structural hamming distance
# In:
# mdl.true  ->  ground truth model (the initial random DAG)
# mdl.est   ->  the estimated model (the structure learning algorithm output)
# Out:
# shd.val,shd.norm   ->   the shd vale and th shd-normalised value
mycomparegraphs3 <- function(mdl.true,mdl.est,mdl.cpdag){
  if(isValidGraph(as(mdl.true,"graphAM")@adjMat,type = "dag")==FALSE)
    warning("'mdl.true' is not a DAG.")
  if(!is.null(mdl.cpdag)){
    #TotalEdges from true model (assumes a DAG)
    amTrue  <- as(mdl.true,"graphAM")@adjMat
    amEst   <- as(getGraph(mdl.est),"graphAM")@adjMat
    amCPDAG <- as(getGraph(mdl.cpdag),"graphAM")@adjMat
    
    uEdges <- get.undiredges(amCPDAG)
    e.list <- uEdges$e.list
    tot.edges <- nrow(e.list)
    #DefinedEdges
    eb.list <- rbind(e.list,e.list[,c(2,1)]) #edgeBidirected.list
    mm<-matrix(FALSE,nrow = nrow(amCPDAG),ncol = ncol(amCPDAG) ) #mask for edges of interest (bidirected)
    mm[eb.list] <- TRUE
    
    amTrue[!(mm)] <- 0
    amCPDAG[!(mm)] <- 0 
    amEst[!(mm)] <- 0 
    mdl.true.am <- graphAM(adjMat = amTrue,edgemode = "directed")
    mdl.est.am <- graphAM(adjMat = amEst,edgemode = "directed")
    
    mdl.true <- as(mdl.true.am,"graphNEL")
    mdl.est  <- as(mdl.est.am,"graphNEL")
  }
  #obtain the number of edges in mdl.true
  total.edges1 <- ncol(edgeMatrix(mdl.true))
  total.edges2 <- sum(as(mdl.true,"graphAM")@adjMat)
  stopifnot(total.edges1 == total.edges2)
  
  shd.val <- shd(mdl.true,mdl.est)
  shd.norm <- shd.val/total.edges1
  if(shd.norm==Inf)
    shd.norm <- 1
  
  return(list(shd.val = shd.val, shd.norm=shd.norm, tot.edg = total.edges1))
}

#==============================================================================
# A function to plot TPR vs TNR
plot.tprfpr <- function(res){
  ss <- subset(res,ALG=="CE-Prop")
  ss_m <- aggregate.data.frame(ss[,1:ncol(ss)-1],by=list(ss$nsamp),FUN = mean)
  plot.new()
  plot(c(0,ss_m[order(ss_m$FPR),]$FPR,1),c(0,ss_m[order(ss_m$TPR),]$TPR,1),type="b",col="blue",xlab="FPR",ylab="TPR",ylim=c(.038,1.0), xlim=c(0.038,1.0),axes=FALSE)
  
  ss <- subset(res,ALG=="PC-")
  ss_m <- aggregate.data.frame(ss[,1:ncol(ss)-1],by=list(ss$nsamp),FUN = mean)
  lines(c(0,ss_m[order(ss_m$FPR),]$FPR,1),c(0,ss_m[order(ss_m$TPR),]$TPR,1),type="b",col="red",xlab="FPR",ylab="TPR",ylim=c(.038,1.0), xlim=c(0.038,1.0),axes=FALSE)
  legend(x =0.65,y=0.4,legend = c("BI+CE","PC"), lty = c(1,1),col=c("blue","red"))
  title(main = "Number of Samples")
  
  axis(1, las=1, at=c(0,0.25,0.5,0.75,1),labels = c("0.0","0.25","0.5","0.75","1.0"))
  axis(2, las=1, at=c(0,0.25,0.5,0.75,1),labels = c("0.0","0.25","0.5","0.75","1.0"))
  
  
  #-----------------------------
  ss <- subset(res,ALG=="CE-Prop")
  ss_m <- aggregate.data.frame(ss[,1:ncol(ss)-1],by=list(ss$p),FUN = mean)
  plot.new()
  plot(c(0,ss_m[order(ss_m$FPR),]$FPR,1),c(0,ss_m[order(ss_m$TPR),]$TPR,1),type="b",col="blue",xlab="FPR",ylab="TPR",ylim=c(.038,1.0), xlim=c(0.038,1.0),axes=FALSE)
  
  ss <- subset(res,ALG=="PC-")
  ss_m <- aggregate.data.frame(ss[,1:ncol(ss)-1],by=list(ss$p),FUN = mean)
  lines(c(0,ss_m[order(ss_m$FPR),]$FPR,1),c(0,ss_m[order(ss_m$TPR),]$TPR,1),type="b",col="red",xlab="FPR",ylab="TPR",ylim=c(.038,1.0), xlim=c(0.038,1.0),axes=FALSE)
  legend(x =0.65,y=0.4,legend = c("BI+CE","PC"), lty = c(1,1),col=c("blue","red"))
  title(main = "Number of Variables")
  
  axis(1, las=1, at=c(0,0.25,0.5,0.75,1),labels = c("0.0","0.25","0.5","0.75","1.0"))
  axis(2, las=1, at=c(0,0.25,0.5,0.75,1),labels = c("0.0","0.25","0.5","0.75","1.0"))
  
  #-----------------------------
  ss <- subset(res,ALG=="CE-Prop")
  ss_m <- aggregate.data.frame(ss[,1:ncol(ss)-1],by=list(ss$BI),FUN = mean)
  plot.new()
  plot(c(0,ss_m[order(ss_m$FPR),]$FPR,1),c(0,ss_m[order(ss_m$TPR),]$TPR,1),type="b",col="blue",xlab="FPR",ylab="TPR",ylim=c(.038,1.0), xlim=c(0.038,1.0),axes=FALSE)
  
  ss <- subset(res,ALG=="PC-")
  ss_m <- aggregate.data.frame(ss[,1:ncol(ss)-1],by=list(ss$BI),FUN = mean)
  lines(c(0,ss_m[order(ss_m$FPR),]$FPR,1),c(0,ss_m[order(ss_m$TPR),]$TPR,1),type="b",col="red",xlab="FPR",ylab="TPR",ylim=c(.038,1.0), xlim=c(0.038,1.0),axes=FALSE)
  legend(x =0.65,y=0.4,legend = c("BI+CE","PC"), lty = c(1,1),col=c("blue","red"))
  title(main = "Background Info")
  
  axis(1, las=1, at=c(0,0.25,0.5,0.75,1),labels = c("0.0","0.25","0.5","0.75","1.0"))
  axis(2, las=1, at=c(0,0.25,0.5,0.75,1),labels = c("0.0","0.25","0.5","0.75","1.0"))
  
}


#==============================================================================
gen.pgm <- function(p){
  # Graphical model
  nodes_name <- paste("V",seq(1,p),sep="")
  dag <- random.graph(nodes = nodes_name, method = "ic-dag", 
        max.in.degree = 2, max.out.degree = 2)
  graphviz.plot(dag)
  # Parametric model (In future change matrix() by list() to allow different number of levels)
  nodes.lv <- matrix(data=c("y","n"),ncol = 2,nrow=p,byrow = TRUE, 
                     dimnames = list(nodes_name,NULL))
  cpt <- list()
  # Instantiated model----- VARIABLE n NO FUNCIONA EN `parents_of_n <- dag$nodes$n$parents`
  for(n in names(dag$nodes)){
    parents_of_n <- dag$nodes[[n]]$parents   # parents of n, e.g. "V1" "V2"
    nparents <- length(parents_of_n)      # eg 2
    parent_nam <- list()                  # to create the list of names for dimnames in the CPT
    for(par in c(n,parents_of_n)){             # to assign the names of levels for each parent
      parent_nam <- append(parent_nam,list(nodes.lv[par,])) # (at the end) eg list(c("y","n"),c("y","n"))
    }
    # to assign the names (V1,V2,...) of each member in parent_nam list, eg list(V1=c("y","n"),V2=c("y","n"))
    names(parent_nam)<-c(n,parents_of_n)
    nlev <- length(nodes.lv[n,])          # number of levels of node n
    
    if(length(parents_of_n)==0){
      dimvector <- c(nlev)
    }else{
      # to generate the (nlev_n,nlev_parents) dimension vector eg. (two parents) c(2,2,2)
      dimvector <- apply(nodes.lv[c(n,parents_of_n),],MARGIN=1,FUN=length) 
    }
    prob <- array(0,dim = as.vector(dimvector),dimnames = parent_nam)
    prob[seq(1,prod(dimvector),by = 2)] <- runif(n=(prod(dimvector)/2))
    prob[seq(2,prod(dimvector),by = 2)] <- 1 - prob[seq(1,prod(dimvector),by = 2)]
    cpt <- append(cpt,list(prob))
  }
  names(cpt) <- nodes_name
  
  return(bn <- custom.fit(dag, cpt))
  
}


#==============================================================================
## Ploting intervals
plot.ace.interv <- function(mat.interv){
  max.intervs <- 5
  df <- data.frame(id =seq(from=1,to=max.intervs,by=max.intervs/nrow(mat.interv)),
                 mid = mat.interv[,3],
                 L = mat.interv[,3]-mat.interv[,4],
                 U = mat.interv[,3]+mat.interv[,4],
                 tail = mat.interv[,1],
                 head = mat.interv[,2])
  p <- ggplot(df, aes(x = mid, y = id)) +
        geom_point(size = 4) +
        geom_errorbarh(aes(xmax = U, xmin = L)) +
        expand_limits(y=c(0,5),x=c(0,max(df$U)+0.5))+
        labs(title = "", x = "ACE", y = "Intervals")+
        scale_y_continuous(breaks=c(1,2,3), 
          labels=c(paste(df$tail[1],"->",df$head[1]),"",
                   paste(df$tail[2],"->",df$head[2])))+
    geom_text(aes(label=paste("[",L," , ",U,"]",sep="")),
              vjust = 0, nudge_y = 0.2)+
    theme(axis.text.x = element_text(colour="black",size=13),
         axis.text.y = element_text(colour="black",size=13),
         axis.title.x = element_text(colour="black",size=14),
         axis.title.y = element_text(colour="black",size=15),
         legend.text = element_text(colour="black",size=12))
  p
  return(p)
}

#==============================================================================
plot.ace.interv2 <- function(mat.interv){
  max.intervs <- 5
  df <- data.frame(id =seq(from=1,to=max.intervs,by=max.intervs/nrow(mat.interv)),
                   mid = mat.interv[,3],
                   L = mat.interv[,3]-mat.interv[,4],
                   U = mat.interv[,3]+mat.interv[,4],
                   tail = mat.interv[,1],
                   head = mat.interv[,2])
  plotCI(y=df$id,x=df$mid,li =df$L ,ui = df$U, err="x", 
         ylim=c(0,max.intervs), xlim=c(-0.1,1.1))
  text(y=(df$id)+0.25,x=(df$mid)-0.1, 
       sprintf("%i->%i: [%.4f - %.4f]",df$tail,df$head,df$L,df$U), 
       cex=1.2, pos=4, col="red") 
  return(NULL)
}

plot.cpdag <- function(am){
  mynet <- network(am,directed = TRUE)
  p <- ggnet2(mynet, node.size = 8, node.color = "gray", 
         edge.size = 0.3, edge.color = "black",
         label = TRUE, arrow.gap = 0.05, arrow.size = 6,
         label.color = "black", mode = "circle")
  return(p)
}

#==============================================================================
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#==============================================================================
orientation <- function(df, beta, cpdag.struct) {
  ## With threshold
  cpdagFlag <- FALSE
  if (beta != -1) {
    Ai.idx <- which(df$Ai >= beta)
    if (length(Ai.idx) > 0) {
      ## sort Ai.idx dscendent to start the orientation according frome the greatest to the lowest Ai value
      Ai.sort <-
        sort(df[Ai.idx, ]$Ai, decreasing = TRUE, index.return = TRUE)
      edgeToadd <-
        matrix(c(df[Ai.sort$ix, ]$y, df[Ai.sort$ix, ]$x),
               nrow = length(Ai.idx),
               byrow = FALSE)
    } else{
      # No edge had an Ai above beta
      edgeToadd <- NULL
    }
  } else{
    ## Without threshold
    Ai.maxidx <- which.max(df$Ai)
    edgeToadd <-
      matrix(c(df[Ai.maxidx, ]$y, df[Ai.maxidx, ]$x), ncol = 2)
  }
  if(!is.null(edgeToadd)){
    cpdag.struct.pdag <- NULL
    #cpdag.struct.pdag <-addBgKnowledge(getGraph(cpdag.struct), edgeToadd[, 1], edgeToadd[, 2])
    if (is.null(cpdag.struct.pdag) == TRUE) {
      #warning("Meek rules could not be added. Trying  adding only one edge.")
      
      cpdag.struct.pdag <-
        removeEdge(as.character(edgeToadd[,2]),
                   as.character(edgeToadd[,1]),
#                   getGraph(cpdag.struct))
                   cpdag.struct)
      
      cpdag.struct.pdag@graphData$edgemode<-"directed"
      
      cpdag.struct.pdag <-
        addEdge(as.character(edgeToadd[,1]),
                as.character(edgeToadd[,2]),
#                getGraph(cpdag.struct.pdag))
                cpdag.struct.pdag)
    }

  }else{
    cpdag.struct.pdag <- cpdag.struct
    cpdagFlag <- TRUE
  }
  
  ## Determine new undefined edges
  undir.res <- get.undiredges(as(cpdag.struct.pdag,"matrix"))
  e.list <- undir.res$e.list
  cpdag.struct <- cpdag.struct.pdag
  
  return(list(cpdag.struct=cpdag.struct, e.list=e.list,cpdagFlag=cpdagFlag))
}


#==============================================================================
real.orient <- function(x,y,myDAG){
  dir.real <- "NaE"
  ceff.real <- NaN
  if(xor(any(myDAG@edgeL[[x]]$edges==y),any(myDAG@edgeL[[y]]$edges==x))==FALSE){
    warning("Problems to determine Real orientation")
  }
  if(any(myDAG@edgeL[[x]]$edges==y)==TRUE){
    dir.real <- sprintf("%i,%i",x,y)
    ceff.real <- causalEffect(myDAG, y, x)
  }else{
    if(any(myDAG@edgeL[[y]]$edges==x)==TRUE){
      dir.real <- sprintf("%i,%i",y,x)
      ceff.real <- causalEffect(myDAG, x, y)
    }
  }
  return(list(dir.real=dir.real, ceff.real=ceff.real))
}