\name{simulate_sde_on_branch}
\alias{simulate_sde_on_branch}
\title{
Simulate a stochastic differential equation on a branch. using
the \pkg{yuima} 
}
\description{
The function simulates a stochastic differential equation on a branch using
the \pkg{yuima} package.
}
\usage{
simulate_sde_on_branch(branch.length, model.yuima, X0, step)
}
\arguments{
  \item{branch.length}{
The length of the branch.
}
  \item{model.yuima}{
A object that \pkg{yuima} can understand in order to simulate
a stochastic differential equation, see Example.
}
  \item{X0}{
The value at the start of the branch.
}
  \item{step}{
The simulation step size that is provided to \pkg{yuima}.
}
}
\details{
The function is a wrapper for calling \code{yuima::simulate()}.

}
\value{
It returns a matrix whose first row are the time points on the 
branch and the remaining rows the values of the trait(s).
}

\references{
Bartoszek, K. and Lio', P (2019). 
Modelling trait dependent speciation with Approximate Bayesian Computation.
Acta Physica Polonica B Proceedings Supplement 12(1):25-47.

Brouste A., Fukasawa M., Hino H., Iacus S. M., Kamatani K., Koike Y., 
Masuda H., Nomura R., Ogihara T., Shimuzu Y., Uchida M., Yoshida N. (2014). 
The YUIMA Project: A Computational Framework for Simulation and Inference of 
Stochastic Differential Equations. Journal of Statistical Software, 57(4): 1-51. 

Iacus S. M., Mercuri L., Rroji E. (2017). 
COGARCH(p,q): Simulation and Inference with the yuima Package. 
Journal of Statistical Software, 80(4): 1-49.
}
\author{Krzysztof Bartoszek}


\seealso{
\code{\link[yuima]{setModel}}, \code{\link[yuima]{setSampling}}, 
\code{\link[yuima]{simulate}},
}

\examples{
## simulate a 3D OUBM process on a branch
set.seed(12345)

A <-c("-(x1-1)-2*x3","-(x2+1)+2*x3",0)
S <- matrix( c( 1, 2, 0, 0, 1 , 0, 0, 0, 
2), 3, 3,byrow=TRUE)
yuima.3d <- yuima::setModel(drift = A, diffusion = S,
state.variable=c("x1","x2","x3"),solve.variable=c("x1","x2","x3") )
X0<-c(0,0,0)
step<-0.5 ## for keeping example's running time short <5s as CRAN policy, 
          ## in reality should be much smaller e.g. step<-0.001          
            
time<-1
simulate_sde_on_branch(time,yuima.3d,X0,step)
}
\keyword{datagen}
