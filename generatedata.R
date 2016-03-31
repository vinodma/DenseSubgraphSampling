library(igraph)
print("hello")
library(igraph)
graph_samps <- seq(2000,8500,500)
sampling_types<-c("RE")
samplingrate<-.22

run_graph = function(iters ){
  iters<-length(graph_samps)
  index <-1
  dt <- data.frame(type="",numNodes="",numEdges="",sampnodes="",sampedges="",Samplingrate="",ActualDensity="",SampledDensity="",ActualCommunitysize="",SampledCommunitysize="",actualmodularity="",sampledmodularity="",Lcommsize_act1="",Lcommsize_samp1 = "",Lcommsize_act2="",Lcommsize_samp2 = "",Lcommsize_act3="",Lcommsize_samp3 = "",Clustcoeff_act="",Clustcoeff_samp="",diameter_act="",diameter_samp="",ks_deg="",pval_deg="",ks_cc="",pval_cc="",ks_pg="",pval_pg="",stringsAsFactors =  F)
  repeat {
    
    flname <- paste0("rmat",graph_samps[index],".txt",sep="")
    
    print(flname)
    dat=read.csv(file(paste0("/home/ubuntu/",flname,sep="")),header = FALSE,sep="\t",skip=3)
    el=as.matrix(dat) # coerc the data into a two-column matrix format that igraph likes
    el[,1]=as.character(el[,1])
    el[,2]=as.character(el[,2])
    g=graph.edgelist(el,directed=FALSE) # turns the edgelist into a 'graph object'
    actedges <- gsize(g)
    actnodes <- vcount(g)
    real_comm <- cluster_walktrap(g)
    real_comm_size <- sizes(real_comm)
    sorted_real_comm_size <- sort(real_comm_size,T)
    Lcommsize_act1<-sorted_real_comm_size[[1]]
    Lcommsize_act2<-sorted_real_comm_size[[2]]
    if(length(sorted_real_comm_size)>2)
      Lcommsize_act3<-sorted_real_comm_size[[3]]
    else
      Lcommsize_act3 <- 0
    Scommsize_act<-min(real_comm_size)
    vlist <- V(g)
    vsize <- length(vlist)
    elist <- as_edgelist(g)
    esize <- nrow(elist)
    wtcactual <- walktrap.community(g)
    actualmodularity = modularity(wtcactual)
    clustcoeff_act=transitivity(g)
    diameter_act <- diameter(g)
    deg_act<-degree(g)
    cc_act<-transitivity(g,type="local")
    pg_act<-page_rank(g)$vector
    g_den <- g
    density_act<-compute_density(g_den)
    
    samplingseq<-seq(.1,.9,.1)
    #samplingrate<-sample(samplingseq,1)
    for(samprate in samplingseq){
      samplingrate <- samprate
    for(i in 1:20){
      samp_type<- sample(sampling_types,1)
      #print(samp_type)
      if(samp_type=="PR")
      {
        g1<-pagerank_sampling(g)
      }
      else if(samp_type=="RE")
      {
        g1<-randomEdgeSampling(g,samplingrate)
      }
      else if(samp_type=="RJ")
      {
        g1<-randomJumpSampling(g)
        
      }
      else
      {
        g1<-randomWalkSampling(g)
      }
      sampnodes<-vcount(g1)
      sampedges<-ecount(g1)
      samplingrate<-sampedges/actedges
      sampled_commg1 <- cluster_walktrap(g1)
      sampled_comm_sizeg1 <- sizes(sampled_commg1)
      sorted_samp_comm_size <- sort(sampled_comm_sizeg1,T)
      
      Lcommsize_samp1<-sorted_samp_comm_size[[1]]
      Lcommsize_samp2<-sorted_samp_comm_size[[2]]
      if(length(sorted_samp_comm_size)>2)
        Lcommsize_samp3<-sorted_samp_comm_size[[3]]
      else
        Lcommsize_samp3 <-0
      Scommsize_samp<-min(sampled_comm_sizeg1)
      
      wtcsampledg1 <- walktrap.community(g1)
      sampledmodularityg1=modularity(wtcsampledg1)
      clustcoeff_samp = transitivity(g1)
      diameter_samp <- diameter(g1)
      deg_samp<-degree(g1)
      cc_samp<-transitivity(g1,type="local")
      pg_samp<-page_rank(g1)$vector
      ksval<-ks.test(deg_act,deg_samp)
      ks_cc<-ks.test(cc_act,cc_samp)
      ks_pg<-ks.test(pg_act,pg_samp)
      density_samp<-compute_density(g1)
      rw <- c(samp_type,actnodes,actedges,sampnodes,sampedges,samplingrate,density_act,density_samp, dim(real_comm_size),dim(sampled_comm_sizeg1)  ,actualmodularity,sampledmodularityg1,Lcommsize_act1,Lcommsize_samp1,Lcommsize_act2,Lcommsize_samp2,Lcommsize_act3,Lcommsize_samp3,clustcoeff_act,clustcoeff_samp,diameter_act,diameter_samp,ksval$statistic,ksval$p.value,ks_cc$statistic,ks_cc$p.value,ks_pg$statistic,ks_pg$p.value)
      dt <- rbind(dt,rw)
      write.csv(dt,file="/home/ubuntu/rmat_pg.csv",append=FALSE)
    }
    }
    iters<- iters-1
    index <- index +1 
    if(iters <1){
      # print(as.data.table(dt))
      dt<-dt[-1,]
      dt$Samplingrate <- as.numeric(dt$Samplingrate)
      return (dt);
    }
    
  }
}


compute_density = function(g){
  gS <- g
  bestaverage <- 0.0
  bestiteration <- 0
  iteration <- 0
  bestgraph <- NULL
  dt<-NULL
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
  return(bestaverage)
}

randomEdgeSampling=function(graph,samplingrate)
{
  #samplingseq<-seq(.1,.9,.1)
  #samplingrate<-sample(samplingseq,1)
  print(samplingrate)
  vlist <- V(graph)
  #size = length(vlist)
  elist <- as_edgelist(graph)
  size=nrow(elist)
  samplededges <- elist[sample(nrow(elist),size=samplingrate*size,replace=FALSE),]
  return(graph.edgelist(samplededges,directed=FALSE) )
  
}

randomJumpSampling<-function(graph,teleport=0.15) {
  
  samplingseq<-seq(.1,1,.2)
  samplingfact<-sample(samplingseq,1)
  ncount<-vcount(graph)*samplingfact
  node <- sample.int(vcount(graph), 1)
  selected <- rep(NA,ncount)
  selected[[1]]<-node
  i<-2 
  
  while(i<=ncount) {
    neigh<-neighbors(graph,node)
    if(length(neigh)==0 | runif(1)<=teleport){
      node <- sample.int(vcount(graph), 1)
    } else {
      node <- sample(neigh,1)
    }
    if (sum(node==selected,na.rm=TRUE)==0) {
      selected[[i]]<-node
      i<-i+1
      #print(paste0("We now have ",i," nodes."))
    }
    
  }
  return(induced.subgraph(graph, selected))
}

pagerank_sampling<-function(g)
{
  samplingseq<-seq(.1,.95,.05)
  samplingfact<-sample(samplingseq,1)
  print(samplingfact)
  vlist <- V(g)
  size = length(vlist)
  #print(samplingfact)
  samplednodes <- sample(V(g),size=samplingfact*size,prob=page_rank(g)$vector,replace=FALSE)
  #print("returned")
  return(induced_subgraph(g,vids=samplednodes)) 
  
}

randomWalkSampling<-function(graph,teleport=0.15) {
  
  samplingseq<-seq(.1,.9,.1)
  samplingfact<-sample(samplingseq,1)
  ncount<-vcount(graph)*samplingfact
  Startnode <- sample.int(vcount(graph), 1)
  selected <- rep(NA,ncount)
  selected[[1]]<-Startnode
  i<-2 
  node<-Startnode
  numiters<-0
  prev_i<--1
  while(i<=ncount) {
    neigh<-neighbors(graph,node)
    # print(length(neigh))
    if(length(neigh)==0 ){
      print ("No neighbors")
    }
    
    if(runif(1)<=teleport){
      node <- Startnode
    } else {
      node <- sample(neigh,1)
    }
    if (sum(node==selected,na.rm=TRUE)==0) {
      selected[[i]]<-node
      
      numiters<-0
      i<-i+1
      prev_i<-i
    }
    else
    {
      numiters<-numiters+1
    }
    print(numiters)
    
    if((i==prev_i)&&(numiters>100))
    {
      Startnode <- sample.int(vcount(graph), 1)
      selected <- rep(NA,ncount)
      selected[[1]]<-Startnode
      i<-2 
      node<-Startnode
      numiters<-0
      prev_i<--1
      print("restarting graph")
    }
    
  }
  return(induced.subgraph(graph, selected))
}
run_graph(0)
