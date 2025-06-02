
library(ggplot2)
library(reshape2)
library(igraph)

load("real_data_final_results/Sequential_computation/Sim_filter_indiv_bands/res_name.Rdata")
n_nodes = dim(G.mat)[1]

#binarize matrix
G.mat = ifelse(G.mat, 1, 0)
G.mat.pop <-  ifelse(G.mat[,1:n_nodes], 1, 0)
G.mat.g <-  ifelse(G.mat[,(n_nodes+1):(n_nodes*2)], 1, 0)

sum(G.mat.pop)
sum(G.mat.g)
#################################

#plot graph population

G_graph.p <- graph_from_adjacency_matrix(G.mat.pop, mode = "undirected", diag = F, weighted=F)
G_graph.p 


layout = layout_in_circle(G_graph.p)
plot(G_graph.p,layout= layout, edge.arrow.size=1, vertex.color="lightgreen", vertex.size=10, 
     
     vertex.frame.color="lightgreen", vertex.label.color="black", 
     
     vertex.label.cex=1, vertex.label.dist=0, edge.curved=0.2, edge.color= "gray", edge.width =E(G_graph.p)$weight, margin=c(0,0,0,0)) 

#plot graph group

G_graph.g <- graph_from_adjacency_matrix(G.mat.g, mode = "undirected", diag = F, weighted=F)
G_graph.g 


layout = layout_in_circle(G_graph.g)
plot(G_graph.g,layout= layout, edge.arrow.size=1, vertex.color="lightgreen", vertex.size=10, 
     
     vertex.frame.color="lightgreen", vertex.label.color="black", 
     
     vertex.label.cex=1, vertex.label.dist=0, edge.curved=0.2, edge.color= "gray", edge.width =E(G_graph.g)$weight, margin=c(0,0,0,0)) 



#### USING OR
G.mat.pop.sim <- G.mat.pop + t(G.mat.pop)
G.mat.pop.sim <- ifelse(G.mat.pop.sim>0, 1, 0)

##### USING AND
G.mat.pop.sim <- G.mat.pop + t(G.mat.pop)
G.mat.pop.sim <- ifelse(G.mat.pop.sim>1, 1, 0)



