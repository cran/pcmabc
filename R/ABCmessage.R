## This file is part of pcmabc

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at bartoszekkj@gmail.com .



.f_printprogress<-function(type,i,dist,par,...){
    switch(type,
	dist=.f_print_dist(i,dist),
	ou1d=.f_print_ou1d(i,dist,par),
	bm1d=.f_print_bm1d(i,dist,par),
    )
    NA
}

.f_print_dist<-function(i,dist,...){print(c(i,dist))}
.f_print_ou1d<-function(i,dist,par,...){cat(paste(i,". dist=",round(dist,digits=7),", a=",round(par$phenotype.model.params$a,digits=7),", s=",round(par$phenotype.model.params$s,digits=7),", s2/2a=",round((par$phenotype.model.params$s^2)/(2*par$phenotype.model.params$a),digits=7),", psi=",round(par$phenotype.model.params$psi,digits=7),"\n",sep=""))}
.f_print_bm1d<-function(i,dist,par,...){cat(paste(i,". dist=",round(dist,digits=7),", s=",round(par$phenotype.model.params$s,digits=7),"\n",sep=""))}

