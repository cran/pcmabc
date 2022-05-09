\name{pcmabc-package}
\alias{pcmabc-package}
\alias{pcmabc}
\docType{package}
\title{
Approximate Bayesian Computations for Phylogenetic Comparative Methods
}
\description{
The package allows for Approximate Bayesian Computations (ABC) inference and simulation 
of stochastic processes evolving on top of a phylogenetic tree. The user is allowed
to define their own stochastic process for the trait(s), be they univariate, multivariate
continuous or discrete. The traits are allowed to influence the speciation
and extinction rates generating the phylogeny. The user provides their own function
that calculates these rates and one of the functions parameters is the current
state of the phenotype. Two functionalities that are missing at the moment
is for the speciation and extinction rates to be time-inhomogenous and 
that speciation events can influence the phenotypic evolution. It is planned to
add this functionality. However, cladogenetic dynamics are possible at the start of the lineage, 
i.e. when a new lineage is separated from the main lineage. Hence, cladogenetic change can only be 
included as an event connected with a lineage (subpopulation) breaking off.


This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
Please understand that there may still be bugs and errors. Use it at your own risk. 
We take no responsibility for any errors or omissions in this package or for any misfortune 
that may befall you or others as a result of its use. Please send comments and report 
bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .
}
\details{
\tabular{ll}{
Package: \tab pcmabc\cr
Type: \tab Package\cr
Version: \tab 1.1.3\cr
Date: \tab 2022-05-06\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
The package allows for Approximate Bayesian Computations (ABC) inference and simulation 
of stochastic processes evolving on top of a phylogenetic tree. The \code{PCM_ABC()}
function is responsible for the inference procedure. The user has to provide their
own functions to simulate the phenotype and tree. The package provides simulation
of the trait under stochastic differential equation (SDE) models by calling \pkg{yuima}.
In this case the user has to provide their own wrapper that creates an object
\pkg{yuima} will be able to handle (see Example). The user also has to provide
a function to calculate the instantaneous birth and death rates (or state that the 
tree is independent of the trait). The package supports some simple birth rate 
functions, see description of \code{PCM_ABC()}.

One is allowed to simulate a trait process and tree process jointly (with the
trait influencing the tree's dynamics). This is done by the function
\code{simulate_phylproc()}. The function \cr \code{simulate_phenotype_on_tree()}
simulates a trait process on top of a provided phylogeny. Finally,
the function \code{simulate_sde_on_branch()} allows one to simulate an SDE
model along a time interval using \pkg{yuima}. 
The trajectory of the trait(s) can be visualized using \code{draw_phylproc()}.
}
\author{
Krzysztof Bartoszek, Pietro Lio'
Maintainer: <krzbar@protonmail.ch>
}
\references{
Bartoszek, K. and Lio', P (2019). 
Modelling trait dependent speciation with Approximate Bayesian Computation.
Acta Physica Polonica B Proceedings Supplement 12(1):25-47.

Kutsukake N., Innan H. (2014) Detecting Phenotypic Selection by Approximate Bayesian 
Computation in Phylogenetic Comparative Methods. In: Garamszegi L. (eds) 
Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. 
Springer, Berlin, Heidelberg

Sheldon R. M. (2006). Simulation. Elsevier Academic Press.
}

\keyword{models}
\keyword{htest}
\keyword{hplot}
\keyword{datagen}


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

sde.params<-list(a11=2.569531,a12=0,a21=0,a22=28.2608,b11=-5.482939,
b21=-34.806936,s11=0.5513215,s12=1.059831,s22=1.247302,s33=1.181376,
psi1=-2.4590422,psi2=-0.6197838)
X0<-c(5.723548,4.103157,3.834698)
step<-0.5 ## for keeping example's running time short <5s as CRAN policy, 
          ## in reality should be much smaller e.g. step<-0.001          
            
simres<-simulate_phenotype_on_tree(phyltree, simulate_mvsl_sde, sde.params, X0, step)

## visualize the simulation
draw_phylproc(simres)

## extract the measurements at the tips
phenotypedata<-get_phylogenetic_sample(simres)

birth.params<-list(scale=1,maxval=2,abcstepsd=0.1,positivevars=c(TRUE,TRUE),
fixed=c(FALSE,TRUE))

sde.params<-list(a11=2.569531,a12=0,a21=0,a22=28.2608,b11=-5.482939,
b21=-34.806936,s11=0.5513215,s12=1.059831,s22=1.247302,s33=1.181376,
psi1=-2.4590422,psi2=-0.6197838,
positivevars=c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE),
abcstepsd=rep(0.1,12))

par0<-list(phenotype.model.params=sde.params,birth.params=birth.params)
fbirth<-"rate_id" ## should be rate_const but used here as an example 
## for birth.params
fdeath<-NULL
X0<-c(5.723548,4.103157,3.834698)
step<-0.5 ## for keeping example's running time short <5s as CRAN policy, 
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
