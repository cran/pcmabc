## This file is part of pcmabc

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .



.f_standard_update<-function(par,par0,baccept,ABCdist,...){
    for (j in 1:length(par)){
	for (k in 1:length(par[[j]])){
	    if ((is.null(par0[[j]]$fixed) || (!par0[[j]]$fixed[k]))){			
		if ((names(par[[j]])[k]!="fixed") && (names(par[[j]])[k]!="abcstepsd") && (names(par[[j]])[k]!="positivevars")){
		    if (!(is.null(par0[[j]]$positivevars)) && (par0[[j]]$positivevars[k]) ) {par[[j]][[k]]<-log(par[[j]][[k]])}
		    if(!baccept){
		        minval<-10;maxval<-10
		        if (!(is.null(par0[[j]]$positivevars)) && (par0[[j]]$positivevars[k]) ) {minval<-log(minval);maxval<-log(maxval)}
		        par[[j]][[k]]<-stats::runif(1,min= (-1)*minval, max=maxval)
		    }
		    else{par[[j]][[k]]<-par[[j]][[k]]+stats::rnorm(1,sd=ABCdist*par0[[j]]$abcstepsd[k])}
		    if (!(is.null(par0[[j]]$positivevars)) && (par0[[j]]$positivevars[k]) ){par[[j]][[k]]<-exp(par[[j]][[k]])}
		}
    	    }
        }
    }
    par
}


.f_unifdraw_update<-function(par,par0,...){
    for (j in 1:length(par)){
	for (k in 1:length(par[[j]])){
	    if ((is.null(par0[[j]]$fixed) || (!par0[[j]]$fixed[k]))){			
		if ((names(par[[j]])[k]!="fixed") && (names(par[[j]])[k]!="abcstepsd") && (names(par[[j]])[k]!="positivevars")){
		    if (!(is.null(par0[[j]]$positivevars)) && (par0[[j]]$positivevars[k]) ) {par[[j]][[k]]<-log(par[[j]][[k]])}
		    minval<-10;maxval<-10
		    if (!(is.null(par0[[j]]$positivevars)) && (par0[[j]]$positivevars[k]) ) {minval<-log(minval);maxval<-log(maxval)}
		    par[[j]][[k]]<-stats::runif(1,min= (-1)*minval, max=maxval)
		    if (!(is.null(par0[[j]]$positivevars)) && (par0[[j]]$positivevars[k]) ){par[[j]][[k]]<-exp(par[[j]][[k]])}
		}
    	    }
        }
    }
    par
}

.f_ou1d_update<-function(par,par0,baccept,ABCdist,phenotypedata,...){
    s2a<-stats::var(phenotypedata[,1])
    adone<-FALSE
    sdone<-FALSE    
    for (j in 1:length(par)){
	for (k in 1:length(par[[j]])){
	    if ((is.null(par0[[j]]$fixed) || (!par0[[j]]$fixed[k]))){			
	        if ((names(par[[j]])[k]!="fixed") && (names(par[[j]])[k]!="abcstepsd") && (names(par[[j]])[k]!="positivevars")){
		    if ((adone || sdone)&&((names(par[[j]])[k]=="a") || (names(par[[j]])[k]=="s"))){
			if (adone){par[[j]][[k]]<-s2a*2*par[[j]]$a+stats::rnorm(1,sd=ABCdist*par0[[j]]$abcstepsd[k])}
			if (sdone){par[[j]][[k]]<-par[[j]]$s/(2*s2a)+stats::rnorm(1,sd=ABCdist*par0[[j]]$abcstepsd[k])}
		    }else{
			if (!(is.null(par0[[j]]$positivevars)) && (par0[[j]]$positivevars[k]) ) {par[[j]][[k]]<-log(par[[j]][[k]])}
	    		    if(!baccept){
				minval<-10;maxval<-10
				if (!(is.null(par0[[j]]$positivevars)) && (par0[[j]]$positivevars[k]) ) {minval<-log(minval);maxval<-log(maxval)}
				par[[j]][[k]]<-stats::runif(1,min= (-1)*minval, max=maxval)
			    }
			    else{par[[j]][[k]]<-par[[j]][[k]]+stats::rnorm(1,sd=ABCdist*par0[[j]]$abcstepsd[k])}
			if (!(is.null(par0[[j]]$positivevars)) && (par0[[j]]$positivevars[k]) ){par[[j]][[k]]<-exp(par[[j]][[k]])}
			if (names(par[[j]])[k]=="a") {adone<-TRUE}
			if (names(par[[j]])[k]=="s") {sdone<-TRUE}
		    }
		}
	    }
	}
    }
    par
}

