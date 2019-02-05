## This file is part of pcmabc

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


PCM_ABC<-function(phyltree,phenotypedata,par0,phenotype.model,fbirth,fdeath=NULL,X0=NULL,step=NULL,abcsteps=200,eps=0.25,fupdate="standard",typeprintprogress="dist",tree.fixed=FALSE,dist_method=c("variancemean","wRFnorm.dist")){
    if (tree.fixed){dist_method[2]<-NA}
    
    if (( is.null(fupdate)) || (is.character(fupdate) && (fupdate=="standard"))){fupdate<-.f_standard_update}
    tree.height<-max(ape::node.depth.edgelength(phyltree))
    if (!is.null(phyltree$root.edge)){tree.height<-tree.height+phyltree$root.edge}
    n<-nrow(phenotypedata) ## we still do not allow for extinct species measurements
    paramsaccepted<-c()
    if (is.null(X0)){## set at average of contemporary
	X0<-apply(phenotypedata,2,mean)
    }
    par<-par0
    for (j in 1:length(par)){
	par[[j]]$fixed<-NULL;par[[j]]$abcstepsd<-NULL;par[[j]]$positivevars<-NULL
    }
    par.estim<-par
    for (j in 1:length(par.estim)){
		for (k in 1:length(par.estim[[j]])){
			    par.estim[[j]][[k]]<-0
		}
	}

    totaldist<-0
    invtotaldist<-0
    
    if (is.character(fbirth)){
	if (fbirth=="rate_id"){fbirth<-.f_rate_id}
	else{if (fbirth=="rate_const"){fbirth<-.f_rate_const}}
    }
    
     if (is.character(fbirth)){
	fbirth<-switch(fbirth,
    	    rate_id=.f_rate_id,
            rate_const=.f_rate_const,
            OU2d_optim_diff=.f_rate_OU2d_optim_diff,
            poly_diff_2=.f_rate_poly_diff_2)
    }
    
    for (i in 1:abcsteps){	
	if (!is.null(par$phenotype.model.params)){simul.params<-par$phenotype.model.params}
	else{stop("No parameters provided to simulate the phenotype (par$phenotype.model.params is NULL).")}
        if (!is.null(par$birth.params)){birth.params<-par$birth.params}
	death.params<-NULL
	if (is.null(fdeath)){n.tips.total<-n}
	else{
	    if (!is.null(par$death.params)){death.params<-par$death.params}
	    n.tips.total<- -1 
	} 
    	
    	ABCdist<-1000+1000*eps
    	baccept<-FALSE
    	tryCatch({
    	    if (!tree.fixed){
    		simdata<-simulate_phylproc(tree.height=tree.height,simul.params=simul.params,X0=X0,fbirth=fbirth,fdeath=fdeath,fbirth.params=birth.params,fdeath.params=death.params,fsimulphenotype=phenotype.model,n.contemporary=n,n.tips.total=n.tips.total,step=step)	
    	    }
    	    if (tree.fixed){
    		simdata<-simulate_phenotype_on_tree(phyltree,fsimulphenotype=phenotype.model,simul.params=simul.params,X0=X0,step=step)
    	    }
    	    simphenotypedata<-get_phylogenetic_sample(simdata)
	    ## tree.fixed not relevant if tree is fixed, then distance is 0 anyway
	    ABCdist<-.calc_ABCdist(simdata$tree,phyltree,simphenotypedata,phenotypedata,dist_method)
	    if (ABCdist < eps){baccept<-TRUE}

	    if (!is.null(typeprintprogress)){
		if (is.function(typeprintprogress)){typeprintprogress(i,ABCdist,par)}
		else{.f_printprogress(typeprintprogress,i,ABCdist,par)}}
	    if(baccept){
		paramsaccepted[[length(paramsaccepted)+1]]<-par
		paramsaccepted[[length(paramsaccepted)]]$distance.from.data<-ABCdist
		totaldist<-totaldist+ABCdist
		if (ABCdist==0){ABCdist<-1e-13}
		paramsaccepted[[length(paramsaccepted)]]$inv.distance.from.data<-1/ABCdist
		invtotaldist<-invtotaldist+1/ABCdist		
	    }
	    par<-fupdate(par,par0,baccept,ABCdist,phenotypedata,phyltree)
	    if (baccept){
		for (j in 1:length(par)){
		    for (k in 1:length(par[[j]])){
			if ((names(par[[j]])[k]!="fixed") && (names(par[[j]])[k]!="abcstepsd") && (names(par[[j]])[k]!="positivevars")){
			    	    par.estim[[j]][[k]]<-par.estim[[j]][[k]]+par[[j]][[k]]*paramsaccepted[[length(paramsaccepted)]]$inv.distance.from.data
			}
		    }
		}
	    }
	},error=function(e){print(paste("Error in simulating phenotype or phylogeny",e))})
    }
    
    for (j in 1:length(par.estim)){
	for (k in 1:length(par.estim[[j]])){
	    if ((is.null(par0[[j]]$fixed) || (!par0[[j]]$fixed[k]))){			
		if ((names(par.estim[[j]])[k]!="fixed") && (names(par.estim[[j]])[k]!="abcstepsd")){
		    par.estim[[j]][[k]]<-par.estim[[j]][[k]]/invtotaldist
		    # http://en.wikipedia.org/wiki/Inverse_distance_weighting
		}
	    }
	    if (!(is.null(par0[[j]]$fixed)) && (par0[[j]]$fixed[k])){			
		if ((names(par.estim[[j]])[k]!="fixed") && (names(par.estim[[j]])[k]!="abcstepsd")){
		    par.estim[[j]][[k]]<-par0[[j]][[k]]
		}
	    }
	}
    }
    list(param.estimate=par.estim, all.accepted=paramsaccepted,sum.dists.from.data=totaldist,sum.inv.dists.from.data=invtotaldist)
}


