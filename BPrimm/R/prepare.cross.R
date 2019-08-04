`prepare.cross` <-
function(y, X, chr, pos=NULL, pheno.name, marker.name, file, sep="\t"){
  if(is.matrix(y)){
    ny = nrow(y)
    n  = ncol(y)
    y  = t(y)
  }else if(is.vector(y)){
    ny = 1
    n  = length(y)
  }else{
    stop("y must be a matrix or a vector")
  }
  
  if(ny != length(pheno.name)){
    stop("dimension of y and pheno.name do not match\n")
  }
  
  if(!is.matrix(X)){
    stop("X must be a matrix\n")
  }
  
  if(nrow(X) != n){
    stop("dimension of X and y do not match\n")
  }
  
  if(ncol(X) != length(marker.name)){
    stop("dimension of X and marker.name do not match\n")
  }
  
  cat("", file=file)
  cat(c(pheno.name, marker.name), file=file, sep=sep, append=TRUE)
  cat("\n", file=file, append=TRUE)
  
  cat(c(rep("", ny), chr), file=file, sep=sep, append=TRUE)
  cat("\n", file=file, append=TRUE)
  
  if(!is.null(pos)){
    cat(c(rep("", ny), pos), file=file, sep=sep, append=TRUE)
    cat("\n", file=file, append=TRUE)
  }

  datf = data.frame(y, X)
  write.table(datf, file = file, append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"))
}

