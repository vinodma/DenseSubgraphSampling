#Calculates dense sub-graph from a given undirected graph

library(igraph)
library(multiplex)

g <- erdos.renyi.game(1000, 1/1000)

gS <- g
bestaverage <- 0.0
bestiteration <- 0
iteration <- 0
bestgraph <- NULL
while(length(V(gS)) > 0){
  grphavg = 2.0*gsize(gS)/length(V(gS))
  
  if(bestaverage <= grphavg){
    bestaverage = grphavg
    bestiteration = iteration
    bestgraph <- gS
  }
  
  degS <- degree(gS, v=V(gS), mode='total')
  dS <- min(degS)
  indexdS <-which.min(degS)
  print(indexdS)
  minS <- V(gS)[indexdS]
  gS <- delete.vertices(gS,minS)
  iteration = iteration + 1
}


print(bestaverage)
