\name{PCM_ABC}
\alias{PCM_ABC}
\title{
ABC estimation for PCMs
}
\description{
The function does parameter estimation through Approximate Bayesian Computations
(ABC) for user defined models of trait and phylogeny evolution. In particular 
the phenotype may influence the branching dynamics.
}
\usage{
PCM_ABC(phyltree, phenotypedata, par0, phenotype.model, fbirth, fdeath = NULL, 
X0 = NULL, step = NULL, abcsteps = 200, eps = 0.25, fupdate = "standard", 
typeprintprogress = "dist", tree.fixed=FALSE, 
dist_method=c("variancemean","wRFnorm.dist"))
}
\arguments{
  \item{phyltree}{
The phylogeny in \code{phylo} format. The tree can be obtained from e.g. a \code{nexus} file 
by the \code{read.nexus()} function from the \pkg{ape} package. The "standard" \pkg{ape} node
indexing is assumed: for a tree with \eqn{n}{n} tips, the tips should have indices \code{1:n} 
and the root index \code{n+1}.
}
  \item{phenotypedata}{
A matrix with the rows corresponding to the tip species while the columns correspond to the traits.
The rows should be in the same as the order in which the 
species are in the phylogeny (i.e. correspond to the node indices \code{1:n}, 
where \eqn{n}{n} is the number of tips). 
}
  \item{par0}{
The starting parameters for the estimation procedure. This is a list of 1, 2 or 3 lists.
The lists have to be named as \cr \code{phenotype.model.params}, \code{birth.params} 
(optional if \code{tree.fixed} is \code{TRUE} and \code{death.params} (optional). 
The \code{phenotype.model.params} list corresponds to parameters governing the 
trait evolution process, the \code{birth.params} to the birth rate of the branching process
and the \code{death.params} to the extinction rate. The last element is optional
as one can have just a pure birth tree. The entries of all the lists should 
be named. Each of the three lists can have three special fields \code{fixed},
\code{abcstepsd}, \code{positivevars}. The field \code{fixed} is a logical
vector of length equalling the number of parameters. If an entry is \code{TRUE},
then that parameter is not optimized over but kept at its initial value
throughout the whole ABC procedure. The field \code{abcstepsd} is a numeric
vector equalling the number of parameters. It is the standard deviation
of the random update of the parameter, i.e. providing control
on how much one wants to jump in each parameter dimension. The field
\code{positivevars} is a logical vector of length equalling the number
of parameters and indicating if the given parameter is to be positive
\code{TRUE} or arbitrary \code{FALSE}. Notice that even if \code{fixed}
has true entries, then corresponding entries have to be present
in \code{abcstepsd} and \code{positivevars} (but their values do not matter).
}
  \item{phenotype.model}{
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
package with the function \code{simulate_sde_on_branch()}. However,
the user needs to write the \code{phenotype.model} function that translates
the vector of parameters into an object that is understandable by \pkg{yuima}
and then call \code{simulate_sde_on_branch()}, see the Example.

The phenotype is simulated prior to the simulation of the speciation/extinction
events on a lineage (see Details). Hence, it is not possible (at the moment) 
to include some special event (e.g. a cladogenetic jump) at branching. Such
dynamics are only possible at the start of the lineage, i.e. when the new lineage 
is separated from the main lineage. Hence, cladogenetic change can only be 
included as an event connected with a lineage (subpopulation) breaking off.
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
possible. The field \code{varnum} indicates the index of the variable
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
  \item{X0}{
The value of the ancestral state at the start of the tree. If \code{NULL}, 
then the mean of the contemporary observations is used.
}
  \item{step}{
The time step size for simulating the phenotype. If not provided,
then calculated as \code{min(0.001,tree.height/1000)}.
}
  \item{abcsteps}{
The number of steps of the ABC search procedure.
}
  \item{eps}{
The acceptance tolerance of the ABC algorithm.
}
  \item{fupdate}{
How should the parameters be updated at each step of the ABC search
algorithm. The user may provide their own function that has to 
handle the following call
\code{fupdate(par,par0,baccept,ABCdist,phenotypedata,phyltree)},
where \code{par} is the current proposal for the parameter,
\code{baccept} a logical variable indicating if \code{par}
were accepted or rejected and \code{ABCdist} the distance
between the observed and simulated under \code{par} data.
The three other parameters \code{par0}, \code{phenotypedata}
and \code{phyltree} are those that were provided in the call
to \code{ABCdist}. The user may write \code{"standard"} in
place of providing a function and then the internal function
\code{.f_standard_update()}. This function makes mean 0 normal jitters
(with standard deviation provided through \code{abcstepsd} from
\code{par0}) for accepted parameters or if they were not accepted
draws new parameters from a uniform on [-10,10] distribution.
}
  \item{typeprintprogress}{
What should be printed out at each step of the ABC search algorithm.
If \code{NULL}, then nothing, otherwise the package offers one possibility,
\code{"dist"}-the iteration number and distance of the simulated dataset from.
The user is free to write their own function here. The first parameter of the
function has to be an integer(the iteration number), the second a real value
(the distance) and the third a list (the proposed model parameters, see \code{par0}
format).
}
  \item{tree.fixed}{
Does the trait value process influence the branching dynamics (\code{FALSE})
or the branching structure is independent of the trait (\code{TRUE}).
}
  \item{dist_method}{
A vector with two entries, the first is the method for calculating   
the distance between the simulated and observed trait data.
The second is is the method for calculating   
the distance between the simulated and observed phylogeny.
Possible values for the phenotype distance are
\code{"variance"}, \code{"variancemean"}, \code{"order"}
and for the distance between phylogenies are
\code{"bdcoeffs"}, \code{"node_heights"}, \code{"logweighted_node_heights"},
\code{"RF.dist"}, \code{"wRF.dist"}, \code{"wRFnorm.dist"}, \code{"KF.dist"},
\code{"path.dist"}, \code{"path.dist.weights"}, \code{"dist.topo.KF1994"},
\cr \code{"treeDist"}, \code{"BHV"} and \code{"BHVedge"}.
See Details for description of the methods. The choice of which distance
calculation method is better seems to depend on the model of evolution,
number of \code{abcsteps} and value of \code{eps}. Some experimentation
is recommended. If the tree is assumed to be fixed, then the tree 
distance method is ignored. 
}
}
\details{
Some details of the function might change. In a future release it should
be possible for the user to provide their own custom distance function,
time-inhomogenous branching and trait simulation. The fields
\code{sum.dists.from.data} and \code{sum.inv.dists.from.data} will 
probably be removed from the output object. 

At the moment the distance function is calculated as
(tree.distance+trait.distance)/2 or only with trait.distance if the tree is assumed 
fixed. The possible distance functions between simulates and observed
phenotypes are
\itemize{
\item{\code{"variance"}}{calculates the Euclidean distance between covariance matrices
estimated from simulated and original data. The differences between entries
are weighted by the sum of their absolute values so that the distance is in [0,1],}
\item{\code{"variancemean"}}{calculates the Euclidean distance between mean vectors and 
covariance matrices estimated from simulated and original data. The differences between entries
are weighted by the sum of their absolute values so that the distance is in [0,1].
The difference between the means and covariances are weighted equally.}
\item{\code{"order"}}{calculates the 
mean squared difference between ordered (by absolute value) tip measurements
(scaled by maximum observed trait in each dimension so that distance is in [0,1]
and then the difference between each pair of traits is scaled by the sum of 
their absolute values to remove effects of scale).}
}
The possible distance functions between simulates and observed
phenotypes are
\itemize{
\item{\code{"bdcoeffs"}}{first using \code{geiger::bd.km()} 
estimates the net diversification rate for both trees and then calculates the distance 
between the trees as the total variation distance between two exponential distributions
with rates equalling the estimated net diversification rates.} 
\item{\code{"node_heights"}}{calculates the root mean square
distance between node heights (scaled by tree height so that distance is in [0,1]
and then the difference between each node height is scaled by the sum of the two heights
to remove effects of scale)}
\item{\code{"logweighted_node_heights"}}{similarly as \code{"node_heights"} but additionally
divides the \cr squared scaled difference is divided by the logarithm of the inverse order (i.e. 
the highest is 1) statistic, to reduce the role of smaller heights. }
\item{\code{"RF.dist"}}{calls \code{phangorn::RF.dist()}}
\item{\code{"wRF.dist"}}{calls \code{phangorn::wRF.dist()} with \code{normalize=FALSE}} 
\item{\code{"wRFnorm.dist"}}{calls \code{phangorn::wRF.dist()} with \code{normalize=TRUE}}
\item{\code{"KF.dist"}}{calls \code{phangorn::KF.dist()}}
\item{\code{"path.dist"}}{calls \code{phangorn::path.dist()} with \code{use.weight=FALSE}}
\item{\code{"path.dist.weights"}}{calls \code{phangorn::path.dist()} with \code{use.weight=TRUE}}
\item{\code{"dist.topo.KF1994"}}{calls \code{ape::dist.topo()} with \code{method="score"}}
\item{\code{"treeDist"}}{calls \code{treespace::treeDist()} with \code{lambda=0}}
\item{\code{"BHV"}}{calls \code{distory::dist.multiPhylo()} with \code{method="geodesic"}}
\item{\code{"BHVedge"}}{calls \code{distory::dist.multiPhylo()} with \code{method="edgeset"}}
}

The tree is simulated by means of a Cox process (i.e. Poisson process
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
\item{param.estimate}{A list, in the form of \code{par0}, with a point
estimate of the parameters. This point estimate is calculated 
from all the accepted points by using inverse distance weighting.
The distances for the weighting are the 
}
\item{all.accepted}{A list with all the accepted parameters during the 
ABC search. This will allow for the presentation of the posterior
distribution of the parameters.}
\item{sum.dists.from.data}{The sum of all the distances between the observed 
data and the simulated data under the accepted parameter sets in the ABC 
search procedure.}
\item{sum.inv.dists.from.data}{The sum of all the inverses of the distances 
between the observed data and the simulated data under the accepted parameter 
sets in the ABC search procedure.}
}
\references{
Sheldon R. M. (2006). Simulation. Elsevier Academic Press.
}
\author{Krzysztof Bartoszek}


\seealso{
\code{\link{simulate_sde_on_branch}}
}
\examples{
## simulate 3d OUBM model
## This example requires the R package TreeSim 
## (to generate the tree in phylo format).
set.seed(12345)

phyltree<-TreeSim::sim.bd.taxa(numbsim=1,n=15,lambda=1,mu=0)[[1]]

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

sde.params<-list(a11=2.569531,a12=0,a21=0,a22=28.2608,b11=-5.482939,
b21=-34.806936,s11=0.5513215,s12=1.059831,s22=1.247302,s33=1.181376,
psi1=10,psi2=-10)
X0<-c(10,4.103157,3.834698)
step<-0.5 ## for keeping example's running time short <5s as CRAN policy, 
            ## in reality should be much smaller e.g. step<-0.001          
simres<-simulate_phenotype_on_tree(phyltree, simulate_mvsl_sde, sde.params, X0, step)

## extract the measurements at the tips
phenotypedata<-get_phylogenetic_sample(simres)

birth.params<-list(rate=10,maxval=10,abcstepsd=0.01,positivevars=c(TRUE,TRUE),
fixed=c(FALSE,TRUE))

sde.params<-list(a11=2.569531,a12=0,a21=0,a22=28.2608,b11=-5.482939,
b21=-34.806936,s11=0.5513215,s12=1.059831,s22=1.247302,s33=1.181376,
psi1=10,psi2=-10,
positivevars=c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE),
abcstepsd=rep(0.1,12))

par0<-list(phenotype.model.params=sde.params,birth.params=birth.params)
fbirth<-"rate_const" 
fdeath<-NULL
X0<-c(10,4.103157,3.834698)
step<-0.05 ## for keeping example's running time short <5s as CRAN policy, 
            ## in reality should be much smaller e.g. step<-0.001          
abcsteps<-2 ## for keeping example's running time short <5s as CRAN policy, 
          ## in reality should be much larger e.g. abcsteps<-500
eps<-1 ## for toy example's output to be useful, 
          ## in reality should be much smaller e.g. eps<-0.25
## estimate parameters
ABCres<-PCM_ABC(phyltree=phyltree,phenotypedata=phenotypedata,
par0=par0,phenotype.model=simulate_mvsl_sde,fbirth=fbirth,fdeath=fdeath,X0=X0,
step=step,abcsteps=abcsteps,eps=eps)
}
\keyword{models}
\keyword{htest}
