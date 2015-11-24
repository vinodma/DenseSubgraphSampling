library(igraph)
library(multiplex)
library(data.table)
set.seed(1234)

samplingrate = 0.8
numnodes = 1000
edgedensity = seq(.4,1,.1)
confidence = .95
epsilon = 2.0
#g <- erdos.renyi.game(1000, 1/1000)

########using file input graph##########
#args <- commandArgs(trailingOnly = TRUE)
#dat=read.csv(file("/Users/vmangipudi/vagrant-xdata/sample.txt"),header = FALSE,sep=" ")
#el=as.matrix(dat) # coerc the data into a two-column matrix format that igraph likes
#el[,1]=as.character(el[,1])
#el[,2]=as.character(el[,2])
#g=graph.edgelist(el,directed=FALSE) # turns the edgelist into a 'graph object'
############################################################

##########using random generated graph###############


compute_density = function(g){
  gS <- g
  bestaverage <- 0.0
  bestiteration <- 0
  iteration <- 0
  bestgraph <- NULL
  while(length(V(gS)) > 0){
    grphavg = 2.0*gsize(gS)/length(V(gS))
    #print(grphavg)
    if(bestaverage <= grphavg){
      bestaverage = grphavg
      bestiteration = iteration
      bestgraph <- gS
    }
    
    degS <- degree(gS, v=V(gS), mode='total')
    dS <- min(degS)
    indexdS <-which.min(degS)
    #print(indexdS)
    minS <- V(gS)[indexdS]
    gS <- delete.vertices(gS,minS)
    iteration = iteration + 1
  }
  
  #print(bestaverage)
  #print(bestiteration)
  return (bestaverage)
}

compute_delta <- function(edges,confidence) {
  value <- (-1.0) *  log(1-confidence) / log(edges)
  return(value)
}

compute_C <- function(nodes, edges, delta,epsilon){
  value <- (12.0  * nodes* (4.0 + delta) * log(edges)) / (epsilon * epsilon)
  return (value)
}

compute_fraction <- function(nodes,edges,confidence,epsilon){
  delta <- compute_delta(edges,confidence)
  return (compute_C( nodes,edges, delta,epsilon)/edges)
}


run_graph = function(iters ){
  index <-0
  idx <- 0
  dt <- data.frame(Iteration="",numNodes="",numEdges="",Samplingrate="",Actualdensity="",Sampleddensity="",ActualCommunitysize="",SampledCommunitysize="",Compare="",stringsAsFactors = F)
  repeat {
    
    #dn <- sample(edgedensity,1)
    dn <- edgedensity[idx+1]
    print(dn)
    g <- erdos.renyi.game(numnodes, dn,directed = F)
    numedges <- gsize(g)
    samplingrate = compute_fraction(numnodes,numedges,confidence,epsilon)
    
    #print(paste0("Sampling rate:",samplingrate ))
    real_density <- compute_density(g)
    real_comm <- cluster_walktrap(g)
    real_comm_size <- sizes(real_comm)
    elist <- as_edgelist(g)
    size = nrow(elist)
    for(i in seq(1,20,1)){
          print (i)
          samplededges <- elist[sample(nrow(elist),size=samplingrate*size,replace=FALSE),]
          g1=graph.edgelist(samplededges,directed=FALSE) 
          sampled_density <- compute_density(g1)
          sampled_comm <- cluster_walktrap(g1)
          sampled_comm_size <- sizes(sampled_comm)
          nodesg <- V(g)
          nodesg1 <- V(g1)
          
          intersect <- nodesg %in% nodesg1
          error_factor <- abs(sampled_density-real_density)/real_density
          #print(error_factor)
          index <- index+1
          report_Str <- paste0("Iteration:", index, " Sampling rate:",samplingrate, " Nodes:",numnodes," Edges:",numedges," Actual density:",real_density," Sampled density:", sampled_density, " Real Community size:", dim(real_comm_size), " Sampled Community size:", dim(sampled_comm_size))
          #print(paste0("Compare: ",compare(real_comm$membership[intersect],sampled_comm$membership,"rand")))
          
          dt <- rbind(dt,c(index,numnodes,numedges,samplingrate, real_density,sampled_density,dim(real_comm_size),dim(sampled_comm_size),compare(real_comm$membership[intersect],sampled_comm$membership,"nmi")))
          
    }
    #print(report_Str)
    iters<- iters-1
    idx <- idx +1
    if(iters <1){
      print(as.data.table(dt))
      dt <- dt[-1,]
      dt$Compare <- as.numeric(dt$Compare)
      dt$Samplingrate <- round(as.numeric(dt$Samplingrate),3)
      dt$Samplingrate <- as.factor(dt$Samplingrate)
      boxplot(dt$Compare ~ dt$Samplingrate,xlab="Sampling Rate",ylab="Compare (community metric)")
      return (dt)
      break;
    }
  }
}
