\name{simulate_phenotype_on_tree}
\alias{simulate_phenotype_on_tree}
\title{
Simulate a stochastic process on top of a phylogeny.
}
\description{
The function simulates phenotypic dataset under user defined models trait 
evolution. The phenotype evolves on top of a fixed phylogeny.
}
\usage{
simulate_phenotype_on_tree(phyltree, fsimulphenotype, simul.params, X0, step)
}

\arguments{
  \item{phyltree}{
The phylogeny in \code{phylo} format. The tree can be obtained from e.g. a \code{nexus} file 
by the \code{read.nexus()} function from the \pkg{ape} package. The "standard" \pkg{ape} node
indexing is assumed: for a tree with \eqn{n}{n} tips, the tips should have indices \code{1:n} 
and the root index \code{n+1}.
}
  \item{fsimulphenotype}{
The name of a function to simulate the phenotype over a period of time. 
The function has to have four parameters (in the following order not by name):
\code{time}, \code{params}, \code{X0} and \code{step}.
The parameter \code{time} is the duration of the simulation, i.e. the length
of the branch. The parameter \code{params} is a list of
the parameters governing the simulation, what will be passed here
is the list \cr \code{phenotype.model.params}, without the fields
\code{fixed}, \code{abstepsd} and \cr \code{positivevars}. It is up
to the function in \code{phenotype.model} to interpret the 
provided parameters. \code{X0} stands for the value at the start of the branch
and \code{step} is a control parameter governing the time 
step size of the simulation. The \pkg{pcmabc} package has inbuilt support for 
simulating the trait as as stochastic differential equation by the \pkg{yuima} 
package with the function \code{simulate_sde_on_branch()}. However
the user needs to write the \code{phenotype.model} function that translates
the vector of parameters into an object that is understandable by \pkg{yuima}
and then call \code{simulate_sde_on_branch()}, see the Example.

The phenotype is simulated prior to the simulation of the speciation/extinction
events on a lineage. Hence, it is not possible (at the moment) 
to include some special event (e.g. a cladogenetic jump) at branching. Such
dynamics are only possible at the start of the lineage, i.e. when the new lineage
is separated from the main lineage. Hence, cladogenetic change can only be 
included as an event connected with a lineage (subpopulation) breaking off.
}
  \item{simul.params}{
The parameters of the stochastic model to simulate the phenotype. 
They should be passed as a named list.
}
  \item{X0}{
The value of the ancestral state at the start of the tree.
}
  \item{step}{
The time step size for simulating the phenotype. If not provided,
then calculated as \code{min(0.001,tree.height/1000)}.
}
}
\value{
\item{tree}{The simulated tree in \code{phylo} format.}
\item{phenotype}{A list of the trajectory of the simulated phenotype. 
Each entry corresponds to a lineage that started at some point of the tree.
Each entry is a matrix with first row equalling time (relative to start of 
the lineage, hence if absolute time since tree's origin is needed it
needs to be recalculated, see \code{draw_phylproc()}'s code) and the next 
rows correspond to the trait value(s).}
\item{root.branch.phenotype}{The simulation on the root branch.
A matrix with first row being time and next rows the simulated trait(s).}
}

\references{
Bartoszek, K. and Lio', P (2019). 
Modelling trait dependent speciation with Approximate Bayesian Computation.
Acta Physica Polonica B Proceedings Supplement 12(1):25-47.
}

\author{Krzysztof Bartoszek}

\seealso{\code{\link{simulate_phylproc}}}

\examples{
## simulate 3d OUBM model
## This example requires the R package ape
## (to generate the tree in phylo format).
set.seed(12345)

phyltree<-ape::rphylo(n=5,birth=1,death=0)
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

birth.params<-list(scale=1,maxval=10000)

sde.params<-list(a11=2.569531,a12=0,a21=0,a22=28.2608,b11=-5.482939,
b21=-34.806936,s11=0.5513215,s12=1.059831,s22=1.247302,s33=1.181376,
psi1=-2.4590422,psi2=-0.6197838)
X0<-c(5.723548,4.103157,3.834698)
step<-0.5 ## for keeping example's running time short <5s as CRAN policy, 
          ## in reality should be much smaller e.g. step<-0.001                      
simres<-simulate_phenotype_on_tree(phyltree, simulate_mvsl_sde, sde.params, X0, step)
}
\keyword{datagen}
