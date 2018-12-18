## This file is part of pcmabc

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at bartoszekkj@gmail.com .



## we assume the tree is binary - no multifurcations at the moment

.create_tree_backbone<-function(tree.height){
## 0 is start of tree
    tree.ape<-vector("list",8)
    names(tree.ape)<-c("edge","tip.label","edge.length","Nnode","node.label" ,"root.edge","node.heights","tree.height")
    tree.ape$edge<-matrix(c(2,1),nrow=1,ncol=2)
    tree.ape$node.heights<-rbind(c(2,tree.height),c(1,0))
    tree.ape$tip.label<-c(1)
    tree.ape$tree.height<-tree.height
    tree.ape$edge.length<-tree.height
    tree.ape$Nnode<-1
    tree.ape$node.label<-c(2)
    tree.ape$root.edge<-0
    class(tree.ape)<-"phylo"
    tree.ape
}

.add_first_split<-function(tree.ape,split.time,branch.length=NULL){
    if ((length(tree.ape$tip.label)!=1)||(tree.ape$edge.length[1]<split.time)){tree.ape<-NA}
    else{
    ## this function is partially generally written, partially specialized for first split (i.e. the branch lengths)
	tree.ape$root.edge<-split.time	
	tree.ape$edge[which(tree.ape$edge>1)]<-tree.ape$edge[which(tree.ape$edge>1)]+1 ## there is only the backbone
	tree.ape$node.heights[which(tree.ape$node.heights[,1]>1),1]<-tree.ape$node.heights[which(tree.ape$node.heights[,1]>1),1]+1
	tree.ape$node.label<-tree.ape$node.label+1
	tree.ape$tip.label<-c(tree.ape$tip.label,max(tree.ape$tip.label)+1)	
	tree.ape$node.heights[1,1]<-0
	tree.ape$node.heights<-	rbind(tree.ape$node.heights,c(max(tree.ape$node.label),tree.ape$edge.length-split.time),c(max(tree.ape$tip.label),0))
	tree.ape$edge<-rbind(tree.ape$edge,c(max(tree.ape$node.label),max(tree.ape$tip.label)))
	tree.ape$edge.length<-rep(tree.ape$edge.length-split.time,2)
	if (!is.null(branch.length)){
	    tree.ape$node.heights[nrow(tree.ape$node.heights),2]<-tree.ape$edge.length[2]-branch.length
	    tree.ape$edge.length[2]<-branch.length
	}
    }
    tree.ape
}

.add_tip_to_branch<-function(tree.ape,branch.number,split.time,branch.length=NULL){
    if((length(tree.ape$edge.length)<branch.number)||(tree.ape$edge.length[branch.number]<split.time)){tree.ape<-NA}
    else{
	numtips<-length(tree.ape$tip.label)    
	tree.ape$edge[which(tree.ape$edge>numtips)]<-tree.ape$edge[which(tree.ape$edge>numtips)]+1 ## there is only the backbone
	tree.ape$node.label<-tree.ape$node.label+1
        tree.ape$node.label<-c(tree.ape$node.label,max(tree.ape$node.label)+1)
	tree.ape$tip.label<-c(tree.ape$tip.label,max(tree.ape$tip.label)+1)	
	tree.ape$node.heights[which(tree.ape$node.heights[,1]>numtips),1]<-tree.ape$node.heights[which(tree.ape$node.heights[,1]>numtips),1]+1
	tree.ape$node.heights<-rbind(tree.ape$node.heights,c(max(tree.ape$node.label),tree.ape$node.heights[which(tree.ape$node.heights[,1]==tree.ape$edge[branch.number,1]),2]-split.time),c(max(tree.ape$tip.label),0))        
	tree.ape$edge<-rbind(tree.ape$edge,c(max(tree.ape$node.label),tree.ape$edge[branch.number,2]),c(max(tree.ape$node.label),max(tree.ape$tip.label)))
        tree.ape$edge[branch.number,2]<-max(tree.ape$node.label)
    	tree.ape$edge.length<-c(tree.ape$edge.length,
    	tree.ape$edge.length[branch.number]-split.time,
	tree.ape$node.heights[nrow(tree.ape$node.heights)-1,2])
	tree.ape$edge.length[branch.number]<-split.time
        if (!is.null(branch.length)){
	    tree.ape$node.heights[nrow(tree.ape$node.heights),2]<-tree.ape$edge.length[length(tree.ape$edge.length)]-branch.length
	    tree.ape$edge.length[length(tree.ape$edge.length)]<-branch.length
	}
	tree.ape$Nnode<-tree.ape$Nnode+1
    }
    tree.ape
}


.mark_bd_events_on_branch<-function(phenotype.trajectory,fbirth,fdeath=NULL,fbirth.params=NULL,fdeath.params=NULL,simul.params=NULL){
    death.rates<-NULL
    death.moment<-NULL
    if (!is.null(fdeath)){
	death.rates<-rbind(phenotype.trajectory[1,],apply(phenotype.trajectory,2,fdeath,params=fdeath.params,sde.params=simul.params)) ## we assume rates not depend on past, only maybe on time
	death.events<-.simulate_nonhomo_Poisson(death.rates)
	if (!is.na(death.events[1])){death.moment<-min(death.events)}
    }
    if (!is.null(death.moment)){phenotype.trajectory<-phenotype.trajectory[,which(phenotype.trajectory[1,]<=death.moment)]}
    birth.rates<-rbind(phenotype.trajectory[1,],apply(phenotype.trajectory,2,fbirth,params=fbirth.params,sde.params=simul.params)) ## we assume rates not depend on past, only maybe on time
    birth.events<-.simulate_nonhomo_Poisson(birth.rates)
    list(births=birth.events,death=death.moment)
}


.simulate_nonhomo_Poisson<-function(rates){
    L<-max(rates[2,])
    events<-c()
    tryCatch({
	Poisson.events<-.simulate_homo_Poisson(max(rates[1,]),L)
	if (length(Poisson.events)>0){
	    event.times<-sapply(Poisson.events,function(x,times,probs){
			p<-probs[max(which(times<=x))]
			sample(c(FALSE,TRUE),1,prob=c(1-p,p))		    
		    },times=rates[1,],probs=rates[2,]/L)
	    events<-Poisson.events[event.times]
	}
    },error=function(e){print(paste("Error in simulating events on branch",e,"rates=",rates))})
    if (length(events)==0){events<-NA}
    events
}

.simulate_homo_Poisson<-function(Tmax,L){
    event.times<-Tmax+1
    tryCatch({
	if (L>0){event.times<-stats::rexp(1,L)}
	s<-5
	while(sum(event.times)<Tmax){
	    event.times<-c(event.times,stats::rexp(as.integer(s*L*Tmax)+1,L))
	    s<-2 
	}
	event.times<-cumsum(event.times)    
    },error=function(e){print(paste("Error in simulating events on branch",e,"Tmax=",Tmax,"L=",L))})
    event.times[which(event.times<=Tmax)]
}

simulate_phylproc<-function(tree.height,simul.params,X0,fbirth,fdeath=NULL,fbirth.params=NULL,fdeath.params=NULL,fsimulphenotype="sde.yuima",n.contemporary=-1,n.tips.total=1000,step=NULL){
## n== -1 so no need to check for NULL
## n.tips.total is a stopping condition

    if (is.character(fbirth)){
	fbirth<-switch(fbirth,
	    rate_id=.f_rate_id,
	    rate_const=.f_rate_const,
	    OU2d_optim_diff=.f_rate_OU2d_optim_diff,
	    poly_diff_2=.f_rate_poly_diff_2)
    }
    
    if (n.contemporary> n.tips.total){n.tips.total = -1}
    if (is.null(step)){step<-min(0.001,tree.height/1000)} 
    if ((!is.function(fsimulphenotype))&&(is.character(fsimulphenotype))&&(fsimulphenotype=="sde.yuima")){fsimulphenotype<-simulate_sde_on_branch}
    root.branch.phenotype<-NULL
    tree.ape<-.create_tree_backbone(tree.height)
    lphenotype<-list()
    lphenotype[[1]]<-.simulate_phenotype_on_branch(tree.height,fsimulphenotype,simul.params,X0,step)
    lspeciations<-list()
    lspeciations[[1]]<-.mark_bd_events_on_branch(lphenotype[[1]],fbirth,fdeath,fbirth.params,fdeath.params,simul.params)$births ## we assume that the backbone cannot die
    next.speciation.branch<-1
    numbranches<-1
    root.branch.phenotype<-c()
    if ((!is.na(lspeciations[[1]][1])) && (n.contemporary!=1) && (n.tips.total!=1)){
	root.branch.phenotype<-lphenotype[[1]][,which(lphenotype[[1]][1,]<=lspeciations[[1]][1])]
        if (is.vector(root.branch.phenotype)){root.branch.phenotype<-matrix(root.branch.phenotype,ncol=1)}
        lphenotype[[1]]<-lphenotype[[1]][,which(lphenotype[[1]][1,]>=lspeciations[[1]][1])]
        if(is.vector(lphenotype[[1]])){lphenotype[[1]]<-matrix(lphenotype[[1]],ncol=1)}
	lphenotype[[1]][1,]<-lphenotype[[1]][1,]-lphenotype[[1]][1,1]
	if(is.vector(lphenotype[[1]])){lphenotype[[1]]<-matrix(lphenotype[[1]],ncol=1)}
	lphenotype[[2]]<-.simulate_phenotype_on_branch(tree.height-lspeciations[[1]][1],fsimulphenotype,simul.params,root.branch.phenotype[-1,ncol(root.branch.phenotype)],step)
	newbds<-.mark_bd_events_on_branch(lphenotype[[2]],fbirth,fdeath,fbirth.params,fdeath.params,simul.params)
	lspeciations[[2]]<-newbds$births
	tree.ape<-.add_first_split(tree.ape,lspeciations[[1]][1],branch.length=newbds$death)
	if((!is.na(lspeciations[[1]][1]))&&(length(lspeciations[[1]])>1)){lspeciations[[1]]<-lspeciations[[1]]-lspeciations[[1]][1];lspeciations[[1]]<-lspeciations[[1]][-1]}
	else{lspeciations[[1]]<-NA}
	
	if ((n.tips.total==2) || ((n.contemporary==2)&&(is.null(newbds$death)))){
	    next.speciation.branch<- 0
	}else{
	    if ((!is.na(lspeciations[[2]][1])) && ( (is.na(lspeciations[[1]][1])) || (lspeciations[[2]][1]<lspeciations[[1]][1])))
		{next.speciation.branch<-2}
	    else{
		if (is.na(lspeciations[[1]][1])){next.speciation.branch<-0}
	    }	
	}

	if (is.null(newbds$death)){num.contemporary.tips<-2}
	else{num.contemporary.tips<-1}
	numbranches<-2
	
	while (next.speciation.branch !=0 ){
	    lphenotype[[numbranches+1]]<-lphenotype[[next.speciation.branch]][,which(lphenotype[[next.speciation.branch]][1,]>=lspeciations[[next.speciation.branch]][1])]
	    if(is.vector(lphenotype[[numbranches+1]])){lphenotype[[numbranches+1]]<-matrix(lphenotype[[numbranches+1]],ncol=1)}
	    lphenotype[[numbranches+1]][1,]<-lphenotype[[numbranches+1]][1,] - lphenotype[[numbranches+1]][1,1]
	    if(is.vector(lphenotype[[numbranches+1]])){lphenotype[[numbranches+1]]<-matrix(lphenotype[[numbranches+1]],ncol=1)}
	    lphenotype[[next.speciation.branch]]<-lphenotype[[next.speciation.branch]][,which(lphenotype[[next.speciation.branch]][1,]<=lspeciations[[next.speciation.branch]][1])]
	    if (is.vector(lphenotype[[next.speciation.branch]])){lphenotype[[next.speciation.branch]]<-matrix(lphenotype[[next.speciation.branch]],ncol=1)}
	    if (ncol(lphenotype[[next.speciation.branch]])>0){
		lphenotype[[numbranches+2]]<-.simulate_phenotype_on_branch(tree.ape$node.heights[which(tree.ape$node.heights[,1]==tree.ape$edge[next.speciation.branch,1]),2]-lspeciations[[next.speciation.branch]][1],fsimulphenotype,simul.params,lphenotype[[next.speciation.branch]][-1,ncol(lphenotype[[next.speciation.branch]])],step)
	    }else{
		j<-next.speciation.branch
		while((ncol(lphenotype[[next.speciation.branch]])==0)&&(j>num.contemporary.tips+1)){
		    j<-which(tree.ape$edge[,2]==tree.ape$edge[j,1])
		}
		if(j>num.contemporary.tips+1){X0<-lphenotype[[j]][-1,ncol(lphenotype[[j]])]}
		else{X0<-root.branch.phenotype[-1,ncol(root.branch.phenotype)]}
		lphenotype[[numbranches+2]]<-.simulate_phenotype_on_branch(tree.ape$node.heights[which(tree.ape$node.heights[,1]==tree.ape$edge[next.speciation.branch,1]),2]-lspeciations[[next.speciation.branch]][1],fsimulphenotype,simul.params,X0,step)
	    }
	    newbds<-.mark_bd_events_on_branch(lphenotype[[numbranches+2]],fbirth,fdeath,fbirth.params,fdeath.params,simul.params)	
	    tree.ape<-.add_tip_to_branch(tree.ape,next.speciation.branch,lspeciations[[next.speciation.branch]][1],branch.length=newbds$death)
	    if((!is.na(lspeciations[[next.speciation.branch]][1]))&&(length(lspeciations[[next.speciation.branch]])>1)){lspeciations[[next.speciation.branch]]<-lspeciations[[next.speciation.branch]]-lspeciations[[next.speciation.branch]][1];lspeciations[[numbranches+1]]<-lspeciations[[next.speciation.branch]][-1]}
	    else{lspeciations[[numbranches+1]]<-NA}
	    lspeciations[[numbranches+2]]<-newbds$births
	    lspeciations[[next.speciation.branch]]<-NA
	    numbranches<-numbranches+2
	    if (is.null(newbds$death)){num.contemporary.tips<-num.contemporary.tips+1}		
	    if ((n.tips.total==(tree.ape$Nnode+1)) || (n.contemporary==num.contemporary.tips)){next.speciation.branch<- 0 }
	    else{	
		next.speciation.branch<- 0
		start.times<-sapply(lspeciations,function(x){res<-NA;if(length(x)>0){res<-x[1]};res})
		if(length(which(is.na(start.times)))<length(start.times)){
		    next.speciation.branch<-which(start.times==min(start.times,na.rm=TRUE))
		}
	    }	
	    
	}
    }else{
	root.branch.phenotype<-lphenotype[[1]]
	lphenotype<-list()
	print("Warning no speciation occurred during simulation")
    }
    list(tree=tree.ape,phenotype=lphenotype,root.branch.phenotype=root.branch.phenotype)
}


