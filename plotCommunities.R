library(sets)

plotCommunities<-function(comm1,comm2,title,vi){
  comm1_size = max(comm1$membership)
  comm2_size = max(comm2$membership)
  matSize = comm1_size + comm2_size
  adjacentMat = matrix(0,matSize,matSize)
  coords = matrix(0,matSize,2)
  
  for(i in 1:min(comm1_size,comm2_size)){
    tmp = as.set(which(comm1$membership %in% i))
    tmp2 = as.set(which(comm2$membership %in% i))
    res = gset_similarity(tmp,tmp2,method = "Jaccard")
    adjacentMat[i,max(comm2_size,comm1_size)+i] = res
    adjacentMat[max(comm2_size,comm1_size)+i,i] = res
  }
    
  g = graph.adjacency(adjacentMat,mode="undirected",weighted = T)
  tmp = max(comm2_size,comm1_size);
  V(g)$type[1:tmp] = TRUE;
  tmp2 = (tmp+min(comm2_size,comm1_size))
  V(g)$type[(tmp+1):tmp2] = FALSE
  pdf(paste(title,".pdf"))
  plot(g, layout=layout.bipartite,vertex.color=c("green","cyan")[V(g)$type+1],edge.width = E(g)$weight * 10)
  title(title)
  text = paste("Communities distance ",round(vi, digits=2));
  mtext(text, side = 1)
  dev.off()
}

