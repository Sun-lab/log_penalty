

# Use the 1st hallmark pathway as a toy example ---
library(GSA)
setwd("/fh/fast/sun_w/yang_liu/mofreds/BRCA")
hp <- GSA.read.gmt("h.all.v6.1.symbols.gmt")
hpNm <- length(hp$geneset.names)
infogene1 <- read.table("data/ens_gene_loci.txt", header = TRUE, 
                        sep = "\t", stringsAsFactors = FALSE)
infogene2 <- read.table("data/Homo_sapiens.GRCh37.66.ensemblID_gene_name.txt", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
genedata <- read.table("data/log_TReC_EA.txt", header = TRUE, 
                       sep = "\t", stringsAsFactors = FALSE)
rownames(genedata) <- genedata[, 1]
genedata <- genedata[, -1]
all(rownames(genedata) %in% infogene2$geneId)
gene_symbols <- infogene2$geneNm[match(rownames(genedata), infogene2$geneId)]
gene_ids <- rownames(genedata)
rownames(genedata) <- gene_symbols
tmp <- hp$genesets[[1]]
tmp2 <- intersect(tmp, rownames(genedata))
variables <- genedata[match(tmp2, rownames(genedata)), ] # observations
variables <- t(variables)


# use the first 100 rows and 100 columns 
variables <- variables[1:100, 1:100]

n <- nrow(variables)
p <- ncol(variables)

rownames(variables) <- seq(1, n, by = 1)
colnames(variables) <- seq(1, p, by = 1)
# end of the data part -----

# load the required packages and functions -----
if (!require("dr")) {
  install.packages("dr") # required for BPrimm
  library(dr)
}  
if (!require("Hmisc")) {
  install.packages("Hmisc") # for hoeffd function in step 2
  library(Hmisc)
}  
if (!require("gpr")){
  install.packages("gpr") # for covNoise function in step 2
  library(gpr)
} 


dir_rp <- "/fh/fast/sun_w/yang_liu/mofreds/Rcodes_2019/"
library("BPrimm", lib.loc = dir_rp)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("/fh/fast/sun_w/yang_liu/mofreds/Rcodes_2019/funs/")

alpha <- 0.01 
recur <- FALSE 

# Step 1 Estimation of the Moral Graph ------

vG_est_BIC <- matrix(0, nrow = p, ncol = p) # BIC
vG_est_EBIC <- matrix(0, nrow = p, ncol = p) # EBIC
vG_est_EBICGG <- matrix(0, nrow = p, ncol = p) # EBICGG

vG_mofreds_BIC <- matrix(0, nrow = p, ncol = p) # BIC
vG_mofreds_EBIC <- matrix(0, nrow = p, ncol = p) # EBIC
vG_mofreds_EBICGG <- matrix(0, nrow = p, ncol = p) # EBICGG

tres <- try({
  
  for(samp_iter in 1:p){
    
    # cat(samp_iter, "=>")
    X <- variables[, -samp_iter]
    y <- variables[, samp_iter]
    
    # compute the max lambda and tau s.t. all variables are 0
    tuning <- gtmGAP(y, X, method = "spline")
    
    pmax <- ncol(X) / 2
    
    # multivariate group-wise adaptive penalization
    # method = "spline"
    # the output for mofreds
    # quadratic spline with one inner knot
    
    I1 <- mGAP3(y, X, method = "spline", 
                lambda = tuning$lambda,
                tau = tuning$tau, pMax = pmax, 
                recursive = recur)
    
    vG_est_BIC[samp_iter, as.numeric(colnames(X)[I1$w_bic])] <- 1 
    vG_est_EBIC[samp_iter, as.numeric(colnames(X)[I1$w_extbic])] <- 1 
    vG_est_EBICGG[samp_iter, as.numeric(colnames(X)[I1$w_extbicgg])] <- 1
    
  }
  
})

if(inherits(tres, "try-error")) next  

vG_mofreds_BIC <- mirror(vG_est_EBIC)
vG_mofreds_EBIC <- mirror(vG_est_EBIC)
vG_mofreds_EBICGG <- mirror(vG_est_EBICGG)

# step 2 Estimation of the skeleton -----
# EBICGG
d <- try({
  temp <- mofredskel(vG_mofreds_EBICGG, variables, VERBOSE= TRUE, 
                     condTest = Zhang_test, dimRed = FALSE, 
                     drmethod = 'sir', trans = 'exp', alpha = alpha)
})

if(inherits(d, "try-error")){
  next
}

fit_exp_EBICGG_zhang_drf_sir <- temp
fit7 <- temp
G_mofreds_EBICGG <- fit7[["G"]]
G_mofreds_EBICGG <- mirror(G_mofreds_EBICGG)

# EBIC
d <- try({
  temp <- mofredskel(vG_mofreds_EBIC, variables, VERBOSE= TRUE, 
                     condTest = Zhang_test, dimRed = FALSE, 
                     drmethod = 'sir', trans = 'exp', alpha = alpha)
})

if(inherits(d, "try-error")){
  next
}

fit_exp_EBIC_zhang_drf_sir <- temp
fit14 <- temp
G_mofreds_EBIC <- fit14[["G"]]
G_mofreds_EBIC <- mirror(G_mofreds_EBIC)

# end ---------

