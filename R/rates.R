## This file is part of pcmabc

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .



.f_rate_const<-function(x,params=NULL,...){
    res<-NA
    if (!is.null(params)){
	if (is.list(params)){
	    if (!is.null(params$rate)){res<-params$rate}
	    if (!is.null(params$switch)){
		if (x[2]>params$switch){res<-params$rate2}
	    }
	}	
	else{
	    if (is.numeric(params)){res<-params}
	}
    }
    if (res<0) {res<-NA}
    res
}

.f_rate_id<-function(x,params=NULL,...){
    res<-x[2]
    if (!is.null(params)){
	if (!is.null(params$varnum)){res<-x[params$varnum]}
	if (!is.null(params$substractbase)){res<-max(0,res-params$substractbase)}
	res<-abs(res)
	if (is.null(params$raise.at.end)||!params$raise.at.end){if (!is.null(params$p)){res<-res^params$p}}
	if (!is.null(params$base)){	    
	    if (!is.null(params$const)){
		if (res<params$base){res<-params$const}
	    }else{if (res<params$base){res<-0}}
	}
	if (!is.null(params$invbase)){	    
	    if (!is.null(params$const)){
		if (res>params$invbase){res<-params$const}
	    }else{if (res>params$invbase){res<-0}}
	}
	if (!is.null(params$scale)){res<-res/params$scale}
	if (!is.null(params$raise.at.end)&&params$raise.at.end){if (!is.null(params$p)){res<-res^params$p}}
	res<-abs(res)
	if (!is.null(params$maxval)){
	    if (is.nan(res)||is.infinite(res)){res<-params$maxval}
	    if(res>params$maxval){res<-params$maxval}
	}
    }
    abs(res)
} ## this just assumes we simulate the rate

.f_rate_poly_diff_2<-function(x,params=NULL,...){
    res<-abs(x[2]-x[3])
    if (!is.null(params)){
	if (!is.null(params$varnum)){res<-abs(x[params$varnum[1]]-x[params$varnum[2]])}
	res<-abs(res)
	if (is.null(params$raise.at.end)||!params$raise.at.end){if (!is.null(params$p)){res<-res^params$p}}
	if (!is.null(params$substractbase)){res<-max(0,res-params$substractbase)}
	if (!is.null(params$base)){	    
	    if (!is.null(params$const)){
		if (res<params$base){res<-params$const}
	    }else{if (res<params$base){res<-0}}
	}
	if (!is.null(params$invbase)){	    
	    if (!is.null(params$const)){
		if (res>params$invbase){res<-params$const}
	    }else{if (res>params$invbase){res<-0}}
	}
	if (!is.null(params$scale)){res<-res/params$scale}
	if (!is.null(params$raise.at.end)&&params$raise.at.end){if (!is.null(params$p)){res<-res^params$p}}
	res<-abs(res)
	if (!is.null(params$maxval)){
	    if (is.nan(res)||is.infinite(res)){res<-params$maxval}
            if(res>params$maxval){res<-params$maxval}
        }
    }
    else{res<-(res)^2}
    abs(res)
} ## this just assumes we simulate the rate


.f_rate_OU2d_optim_diff<-function(x,params=NULL,sde.params=NULL){
    res<-0
    if (!is.null(params)&&!is.null(sde.params)	){
	if (!is.null(sde.params)){
	    if (is.null(params$revpredictor) || !params$revpredictor){
		if (!(is.null(sde.params$a11)||is.null(sde.params$a12)||is.null(sde.params$psi1)||is.null(sde.params$psi2))){
		    res<-abs(x[2]-((-1)*sde.params$a12*x[3]/sde.params$a11+sde.params$psi1+sde.params$a12*sde.params$psi2/sde.params$a11))
		}
	    }else{
		if (!(is.null(sde.params$a22)||is.null(sde.params$a21)||is.null(sde.params$psi1)||is.null(sde.params$psi2))){
		    res<-abs(x[3]-((-1)*sde.params$a21*x[2]/sde.params$a22+sde.params$psi2+sde.params$a21*sde.params$psi1/sde.params$a22))
		}
	    }	    
	}
	if (!is.null(params$substractbase)){res<-max(0,res-params$substractbase)}
	if (!is.null(params$base)){	    
	    if (!is.null(params$const)){
		if (res<params$base){res<-params$const}
	    }else{if (res<params$base){res<-0}}
	}
	if (!is.null(params$invbase)){	    
	    if (!is.null(params$const)){
		if (res>params$invbase){res<-params$const}
	    }else{if (res>params$invbase){res<-0}}
	}
	if (!is.null(params$scale)){res<-res/params$scale}
	if (!is.null(params$p)){res<-res^params$p}else{res<-res^2}
	res<-abs(res)
	if (!is.null(params$maxval)){
	    if (is.nan(res)||is.infinite(res)){res<-params$maxval}
            if(res>params$maxval){res<-params$maxval}
        }
    }
    abs(res)
} ## this just assumes we simulate the rate


