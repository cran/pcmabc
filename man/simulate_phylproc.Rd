\name{simulate_phylproc}
\alias{simulate_phylproc}
\title{
Simultaneously simulate a phylogeny and a trait evolving on it.
}
\description{
The function does simulates a phylogeny and phenotypic dataset under 
user defined models of trait and phylogeny evolution. In particular 
the phenotype may influence the branching dynamics.
}
\usage{
simulate_phylproc(tree.height, simul.params, X0, fbirth, fdeath=NULL, 
fbirth.params=NULL, fdeath.params=NULL, fsimulphenotype="sde.yuima", 
n.contemporary=-1, n.tips.total=1000, step=NULL)
}
\arguments{
  \item{tree.height}{
The height of the desired output tree. The simulated tree 
is conditioned to be of a certain height. 
}
  \item{simul.params}{
The parameters of the stochastic model to simulate the phenotype. 
They should be passed as a named list. 
}
  \item{X0}{
The value of the ancestral state at the start of the tree.
}
  \item{fbirth}{
A function that returns the birth rate at a given moment.
The fist parameter of the function has to correspond to the value of the
phenotype (it is a vector first element is time and the others the values of the trait(s)) 
and the second to the list of birth parameters \code{birth.params},
see \code{par0}. The time entry of the phenotype vector is at the moment relative to 
an unspecified (from \code{fbirth}'s perspective) speciation event on the phylogeny.
Hence, it cannot be used as for writing a time-inhomogenous speciation
function. The speciation process is assumed to be time homogeneous in the 
current implementation. The package has support for two inbuilt rate functions.
The can be indicated by passing a string in \code{fbirth}: either \code{"rate_id"}
or \code{"rate_const"}. The string \code{"rate_const"} corresponds to a 
constant rate and has to have the rate's value in field called \code{$rate}.
However, a switching of rates is allowed. If the value of the first trait
exceeds a certain threshold (provided in field \code{$switch} of 
birth parameters), then the rate is changed to the value in \code{$rate2},
see body of hidden function \code{.f_rate_const()}, in file \code{rates.R}.
The string \code{"rate_id"} corresponds to the \code{.f_rate_id()} function, 
in file \code{rates.R}. If the birth parameters are \code{NULL}, then the
rate equals the value of the first phenotype. However, a number 
of linear, threshold and power transformations of the rate are
possible.The field \code{varnum} indicates the index of the variable
to take as the one influencing the rate (remember to add 1 for the time entry).
Then, if \eqn{x}{x} stands for the trait influencing the branching rate it is 
transformed into a rate by the following fields in the following order. 
Set \code{rate<-x} and let \code{params} correspond to the list containing
the branching parameters.
\itemize{
\item{\code{substractbase}}{rate<-max(0,params$rate-substractbase)}
\code{rate<-abs(rate)}
\item{\code{p} and if \code{is.null(params$raise.at.end) || !params$raise.at.end}}{
rate<-rate^params$p}
\item{\code{base} and \code{!is.null(params$const)}} \cr {\code{if (res<params$base){rate<-params$const}}}
\item{\code{base} and \code{is.null(params$const)}} \cr {\code{if (res<params$base){rate<-0}}}
\item{\code{invbase} and \code{!is.null(params$const)}} \cr {\code{if (res>params$invbase){rate<-params$const}}}
\item{\code{invbase} and \code{is.null(params$const)}} \cr {\code{if (res>params$invbase){rate<-0}}}
\item{\code{scale}}{rate<-rate/params$scale}
\item{\code{p} and if \code{params$raise.at.end}}{rate<-rate^params$p}
\code{res<-abs(res)}
\item{\code{maxval}}{\code{if(rate>params$maxval){rate<-params$maxval}}}
}
}
  \item{fdeath}{
A function that returns the birth rate at a given moment. Its structure
is the same as \code{fbirth}. The current version of the  package does 
not provide any support for any inbuilt function.
}
  \item{fbirth.params}{
The parameters of the birth rate, to be provided to \code{fbirth}.
They have to be a named list. 
}
  \item{fdeath.params}{
The parameters of the death rate, to be provided to \code{fdeath}.
They have to be a named list. 
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
events on a lineage (see Details). Hence, it is not possible (at the moment) 
to include some special event (e.g. a cladogenetic jump) at branching. Such
dynamics are only possible at the start of the lineage, i.e. when the new lineage
separated from the main lineage. Hence, cladogenetic change can only be 
included as an event connected with a lineage (subpopulation) breaking off.
}
  \item{n.contemporary}{
The number of contemporary species to generate. If equals -1, then ignored.
Otherwise when the tree reaches \code{n.contemporary} tips at height
\code{tree.height} the simulation is stopped. However, there is no
conditioning on this value, it is just an upper bound. It may
just happen that due to the birth and death rate functions the
process stops before reaching the target number of tips.
}
  \item{n.tips.total}{
The total (contemporary and extinct) number of tips to generate.
If equals -1, then ignored. Otherwise when the tree reaches 
\code{n.tips.total} tips the simulation is stopped. However, there is 
no conditioning on this value, it is just an upper bound. It may
just happen that due to the birth and death rate functions the
process stops before reaching the target number of tips.
}
  \item{step}{
The time step size for simulating the phenotype. If not provided,
then calculated as \code{min(0.001,tree.height/1000)}.
}
}
\details{
The tree is simulate by means of a Cox process (i.e. Poisson process
with random rate). First the trait is simulated along the spine of a 
tree, i.e. a lineage of duration \code{tree.height}. Then, along this
spine the birth and death rates are calculated (they may depend 
on the value of the phenotype). The maximum for each rate is calculated
and a homogeneous Poisson process with the maximum rate is simulated.
Then, these events are thinned. Each event is retained with probability
equalling true rate divided by maximum of rate (p. 32, Sheldon 2006). 
All speciation events are retained until the first death event.
}
\value{
\item{tree}{The simulated tree in \code{phylo} format.}
\item{phenotype}{A list of the trajectory of the simulated phenotype. 
The i-th entry of the list corresponds to the trait's evolution on
the i-th edge (as in i-th row of \code{phyltree$edge}) of the tree.
Each entry is a matrix with first row equalling time (relative to the 
start of the branch) and the next rows correspond to the trait value(s).}
\item{root.branch.phenotype}{The simulation on the root branch.
A matrix with first row being time and next rows the simulated trait(s).}
}
\references{
Bartoszek, K. and Lio', P (2018). 
Modelling trait dependent speciation with Approximate Bayesian Computation.
arXiv:1812.03715

Sheldon R. M. (2006). Simulation. Elsevier Academic Press.
}
\author{Krzysztof Bartoszek}


\seealso{\code{\link{PCM_ABC}}}

\examples{
## simulate 3d OUBM model with id branching rate
set.seed(12345)

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

birth.params<-list(scale=1,maxval=2)

sde.params<-list(a11=2.569531,a12=0,a21=0,a22=28.2608,b11=-5.482939,
b21=-34.806936,s11=0.5513215,s12=1.059831,s22=1.247302,s33=1.181376,
psi1=-2.4590422,psi2=-0.6197838)
X0<-c(5.723548,4.103157,3.834698)
step<-0.25 ## for keeping example's running time short <5s as CRAN policy, 
            ## in reality should be much smaller e.g. step<-0.001          
            
simres<-simulate_phylproc(3.5, sde.params, X0, fbirth="rate_id", fdeath=NULL, 
fbirth.params=NULL, fdeath.params=NULL, fsimulphenotype=simulate_mvsl_sde, 
n.contemporary=5, n.tips.total=-1, step=step)
}
\keyword{datagen}
