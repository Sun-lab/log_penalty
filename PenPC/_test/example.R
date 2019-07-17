library(PEN, lib.loc="/Users/minjinha/Rlibs/") 

# ================================= #
#           connectedComp           #
# ================================= #

# simulate a DAG following ER model
G = simul.ER(p=100,pE=0.02,n=30)$G
Gmat = as(G,"matrix")
Gmat[Gmat!=0|t(Gmat)!=0] = 1 # skeleton of the DAG G
(1:nrow(Gmat))[connectedComp(G=Gmat,x=c(1,2))] # indices for connected components of c(1,2)

# ================================= #
#           ne.PEN                  #
# ================================= #


# simulate a DAG following ER model
p=100
simul = simul.ER(p=p,pE=0.02,n=30)
dat =  simul$X
# Estimate neighbors of all vertices
coeff = ne.PEN(dat=dat,nlambda=100,ntau=10,V=1:p,order=TRUE,verbose=TRUE)
# Fit Partial correlation graph with OR rule
Gmat = matrix(0,p,p)
Gmat[coeff!=0 | t(coeff)!=0] =1
PCgraph = as(Gmat,"graphNEL")


# ================================= #
#           simul.BA                 #
# ================================= #
p = 100 # number of vertices
e = 1 # one edge added in each time step
n = 30 # sample size
simul=simul.BA(p,e,n)


# ================================= #
#           simul.ER                 #
# ================================= #
p = 100 # number of vertices
pE = 0.02 # edge inclusion probability
n = 30 # sample size
simul=simul.ER(p,pE,n)


# ================================= #
#           skeleton.stable         #
# ================================= #
alpha = 0.01
p = 100 # number of vertices
e = 1 # one edge added in each time step
n = 30 # sample size
simul=simul.BA(p,e,n)

dat = simul$X
indepTest = gaussCItest
suffStat  = list(C = cor(dat), n = n) 

fit.pc.stable = skeleton.stable(suffStat, indepTest, p, alpha)



# ================================= #
#           skeletonPENstable       #
# ================================= #
alpha = 0.01 # significance level for a partial correlation testing
p = 100 # number of vertices
e = 1 # one edge added in each time step
n = 30 # sample size
simul=simul.BA(p,e,n)

dat = simul$X
indepTest = gaussCItest
suffStat  = list(C = cor(dat), n = n) 

coeff = ne.PEN(dat=dat,nlambda=100,ntau=10,V=1:p,order=TRUE,verbose=TRUE)
edgeWeights = matrix(0,p,p)
edgeWeights[coeff!=0|t(coeff)!=0] =1

fit.penpc  = skeletonPENstable(suffStat, indepTest, as.integer(p), alpha, 
edgeWeights=edgeWeights, verbose=TRUE)

