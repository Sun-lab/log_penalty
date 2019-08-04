`print.MAL` <-
function(x, ...)
{

  ww = which(names(x) != "yt")
  
  for(i in ww){
    cat(names(x)[i], ":\n", sep="")
    print(x[[i]])
    cat("\n")
  }
  
}

