`mirror` <- function(graph){
  graph2 <- graph
  for (i in 1:nrow(graph)){
    for (j in 1:nrow(graph)){
      if(graph[i,j] == 1){
        graph2[j,i] <- 1
      }
    }
  }
  
  return(graph2)
}
