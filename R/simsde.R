## This file is part of pcmabc

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .



simulate_sde_on_branch<-function(branch.length,model.yuima,X0,step){
    numsteps<-branch.length/step
    options(warn= -1) ## we do not want warning output of yuima 
    samp<-yuima::setSampling(Terminal=branch.length, n=numsteps)
    simulobj.yuima <- yuima::setYuima(model=model.yuima, sampling=samp)
    simulobj.yuima <- yuima::simulate(simulobj.yuima, xinit=X0,space.discretized=TRUE)
    time.grid.length<-length(simulobj.yuima@sampling@grid[[1]])
    options(warn= 0)
    
    sdedim<-length(X0)
    
    trait_data<-NA
    if (is.element(".Data",methods::slotNames(simulobj.yuima@data@original.data))){
	trait_data<-simulobj.yuima@data@original.data@.Data
	time_data<-simulobj.yuima@sampling@grid[[1]]
    }else{
	if (length(X0)==1){
	    time_data<-attributes(simulobj.yuima@data@original.data)$index ## get at indices of zoo object
	    trait_data<-c(simulobj.yuima@data@original.data)
	    class(trait_data)<-"numeric"
	    names(trait_data)<-NULL
	    
	    vtorem<-which(time_data>branch.length)
	    if (length(vtorem)>0){
		trait_data<-trait_data[-vtorem]
		time_data<-time_data[-vtorem]
	    }
	    time.grid.length<-length(time_data)
	}else{
	    stop("Cannot extract simulated by yuima phenotypic data on branch")
	}	
    }
    timepoints<-length(trait_data)/sdedim
    if (timepoints!=time.grid.length){warning("Error in Yuima, wrong lengths of grid")}
    rbind(time_data[1:timepoints],matrix(trait_data,ncol=timepoints,byrow=TRUE))
}

.simulate_phenotype_on_branch<-function(branch.length,fsimulphenotype,simul.params,X0,step){
    res<-matrix(c(0,X0),ncol=1)
    if (branch.length>=step){
	systime<-as.numeric(Sys.time()); set.seed((systime - floor(systime)) * 1e8); 
	res<-fsimulphenotype(branch.length,simul.params,X0,step)
    }
    res
}


.f_phenotype_variance_on_branch<-function(branch.length,simul.params,X0,step){
    X0<-X0[1]
    mphenotype<-simul.params$fsimulphenotype(branch.length,simul.params$funcparams,X0,step)
    if (simul.params$n>1){
	lphenotype<-replicate(simul.params$n-1,simul.params$fsimulphenotype(branch.length,simul.params$funcparams,X0,step)[-1,],simplify=FALSE)
	for (i in 1:length(lphenotype)){
	    mphenotype<-rbind(mphenotype,lphenotype[[1]])
	}
    }
    rbind(mphenotype[1,],apply(mphenotype[-1,],2,mean),apply(mphenotype[-1,],2,var))
}

simulate_phenotype_on_tree<-function(phyltree,fsimulphenotype,simul.params,X0,step){
    if ((!is.function(fsimulphenotype))&&(is.character(fsimulphenotype))&&(fsimulphenotype=="sde.yuima")){fsimulphenotype<-simulate_sde_on_branch}
    numbranch<-nrow(phyltree$edge)
    root.branch.phenotype<-matrix(c(0,X0),ncol=1)
    if (!is.null(phyltree$root.edge)){
	root.branch.phenotype<-.simulate_phenotype_on_branch(phyltree$root.edge,fsimulphenotype,simul.params,X0,step)
	X0<-root.branch.phenotype[-1,ncol(root.branch.phenotype)]
    }else{
	phyltree$root.edge<-0
    }
    lphenotype=vector("list",numbranch)
    vBranchToSim<-rep(1,numbranch) ## we do not assume that we can have a forest
    RootId<-length(phyltree$tip.label)+1
    vEdgesFromRoot<-which(phyltree$edge[,1]==RootId)
    mBranchQueue<-matrix(c(vEdgesFromRoot[1],X0),ncol=1,nrow=1+length(X0))
    if (length(vEdgesFromRoot)>1){
	for (j in 2:length(vEdgesFromRoot)){
	    mBranchQueue<-cbind(mBranchQueue,c(vEdgesFromRoot[j],X0))
	}
    }
    while(sum(vBranchToSim)>0){
	EdgeToSim<-mBranchQueue[1,1]
	lphenotype[[EdgeToSim]]<-.simulate_phenotype_on_branch(phyltree$edge.length[EdgeToSim],fsimulphenotype,simul.params,mBranchQueue[-1,1],step)
	vEdgesFromBranch<-which(phyltree$edge[,1]==phyltree$edge[EdgeToSim,2])
	if (length(vEdgesFromBranch)>0){
	    for (Edge in vEdgesFromBranch){
		mBranchQueue<-cbind(mBranchQueue,c(Edge,lphenotype[[EdgeToSim]][-1,ncol(lphenotype[[EdgeToSim]])]))
	    }
	}
	if (ncol(mBranchQueue)==1){mBranchQueue<-NULL}
	else{
	    bOneBranch<-FALSE
	    if (ncol(mBranchQueue)==2){bOneBranch<-TRUE}
	    mBranchQueue<-mBranchQueue[,-1]
	    if (bOneBranch){mBranchQueue<-matrix(mBranchQueue,ncol=1)}
	}
	vBranchToSim[EdgeToSim]<-0
    }
    ## "Correct" ape tree too include the additional fields that we require in PhyloSDE
    bHeightRev<-FALSE
    if (is.null(phyltree$node.heights)){
	phyltree$node.heights<-cbind(1:(phyltree$Nnode+length(phyltree$tip.label)),node.depth.edgelength(phyltree)+phyltree$root.edge)
	bHeightRev<-TRUE
    }
    if (is.null(phyltree$tree.height)){phyltree$tree.height<-phyltree$node.heights[1,2]}
    if (bHeightRev){phyltree$node.heights[,2]<-phyltree$tree.height-phyltree$node.heights[,2]}
    if (is.null(phyltree$node.label)){phyltree$node.label<-RootId:(RootId+phyltree$Nnode)}
    list(tree=phyltree,phenotype=lphenotype,root.branch.phenotype=root.branch.phenotype)
}
