\name{get_phylogenetic_sample}
\alias{get_phylogenetic_sample}
\title{
Retrieve contemporary sample
}
\description{
The function retrieved the contemporary sample from an object
simulated by \pkg{pcmabc}.
}
\usage{
get_phylogenetic_sample(pcmabc_simulobj, bOnlyContemporary=FALSE, 
tol=.Machine$double.eps^0.5)
}
\arguments{
  \item{pcmabc_simulobj}{
The output simulated object by \code{simulate_phylproc()}
\cr or \code{simulate_phenotype_on_tree()}.
}
  \item{bOnlyContemporary}{
Logical, should all tip measurements be extracted (\code{FALSE})
or only those corresponding to non-extinct tips (\code{TRUE}).
In the latter case be careful with the output (see Value).
}
  \item{tol}{
The tolerance to check if the tip is contemporary, i.e.
the tip is distant from the root by tree.height +/- \code{tol}.
}
}
\value{
The function returns a matrix with rows corresponding to tips.
If there are extinct species, but \code{bOnlyContemporary} was
set to \code{TRUE}, then the extinct lineages have to be removed
from the tree before any further analysis. Also after removing
extinct lineages one has to make sure that the order of the 
contemporary nodes did not change in the tree. Otherwise, 
the rows of the output matrix will not correspond to the
indices of the tree's tips.
}

\references{
Bartoszek, K. and Lio', P (2018). 
Modelling trait dependent speciation with Approximate Bayesian Computation.
arXiv:1812.03715
}

\author{Krzysztof Bartoszek}



\examples{
## simulate 3d OUBM model
## This example requires the R package TreeSim 
## (to generate the tree in phylo format).
set.seed(12345)

phyltree<-TreeSim::sim.bd.taxa(numbsim=1,n=5,lambda=1,mu=0)[[1]]
simulate_mvsl_sde<-function(time,params,X0,step){
    A <- c(paste("(-",params$a11,")*(x1-(",params$psi1,"))
    -(",params$a12,")*(x2-(",params$psi2,"))-(",params$b11,")*x3",sep=""),
    paste("(-",params$a21,")*(x1-(",params$psi1,"))
    -(",params$a22,")*(x2-(",params$psi2,"))-(",params$b21,")*x3",sep=""),0)
    S <- matrix( c( params$s11, params$s12, 0, 0, params$s22 
    , 0, 0, 0, params$s33), 3, 3,byrow=TRUE)
    yuima.3d <- yuima::setModel(drift = A, diffusion = S,
    state.variable=c("x1","x2","x3"),solve.variable=c("x1","x2","x3") )
    simulate_sde_on_branch(time,yuima.3d,X0,step)
}
phenotype.model<-simulate_mvsl_sde

birth.params<-list(scale=1,maxval=10000,abcstepsd=1,positivevars=c(TRUE),
fixed=c(FALSE,TRUE,TRUE,TRUE,TRUE))

sde.params<-list(a11=2.569531,a12=0,a21=0,a22=28.2608,b11=-5.482939,
b21=-34.806936,s11=0.5513215,s12=1.059831,s22=1.247302,s33=1.181376,
psi1=-2.4590422,psi2=-0.6197838)
X0<-c(5.723548,4.103157,3.834698)
step<-0.5 ## for keeping example's running time short <5s as CRAN policy, 
          ## in reality should be much smaller e.g. step<-0.001          
simres<-simulate_phenotype_on_tree(phyltree, simulate_mvsl_sde, sde.params, X0, step)
mPhylSamp<-get_phylogenetic_sample(simres)
}
\keyword{datagen}
