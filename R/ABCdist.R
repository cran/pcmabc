## This file is part of pcmabc

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .



.calc_phyl_dist<-function(tree1,tree2,method="node_heights"){
## tree1 simtree, tree2 truetree
    tree_dist<-0
    
    if (!is.na(method)){
	if (method!="bdcoeffs"){
	    tryCatch({tree1<-ape::drop.fossil(tree1)},error=function(e){e<-NULL;e})
	    tryCatch({tree2<-ape::drop.fossil(tree2)},error=function(e){e<-NULL;e})
        }
        
        if ((method!="bdcoeffs") && (method!="node_heights")){
    	    ## the other methods require the trees to have the same number of tips
    	    ## as the tips are added here sequentially, then it is the ones with
    	    ## the earlier labels that seem more important in the generated sample
    	    common_tip_labels<-intersect(1:length(tree1$tip.label),1:length(tree2$tip.label))
	    
	    vtips_del1<-setdiff(1:length(tree1$tip.label),common_tip_labels)
    	    if (length(vtips_del1)>0){tree1<-ape::drop.tip(tree1,vtips_del1)}
    	    
    	    vtips_del2<-setdiff(1:length(tree2$tip.label),common_tip_labels)
    	    if (length(vtips_del2)>0){tree2<-ape::drop.tip(tree2,vtips_del2)}
    	    if (!is.null(tree1) && !is.null(tree2)) {tree1$tip.label<-tree2$tip.label}
    	}
	if (is.null(tree1) || is.null(tree2)) {return(NA)}
	
	if (method=="BHV"){
	    if (requireNamespace("distory",quietly=TRUE)){ 
		ltre<-list(tree1,tree2)
		class(ltre)<-"multiPhylo"
		tree_dist<-distory::dist.multiPhylo(ltre,method="geodesic")[[1]] 
	    }else{method<-"wRFnorm.dist"}
	}
	if (method=="BHVedge"){
	    if (requireNamespace("distory",quietly=TRUE)){ 
		ltre<-list(tree1,tree2)
		class(ltre)<-"multiPhylo"
		tree_dist<-distory::dist.multiPhylo(ltre,method="edgeset")[[1]] 
	    }else{method<-"RF.dist"}
	}    

	if ((method=="node_heights")||(method=="logweighted_node_heights")){
	    s1<-sort(ape::node.depth.edgelength(tree1),decreasing = TRUE)
	    s2<-sort(ape::node.depth.edgelength(tree2),decreasing = TRUE)
	    s1<-s1/s1[1]
	    s2<-s2/s2[1]
	    n<-min(length(s1),length(s2))

	    vtokeep<-which(sapply(s1[1:n]+s2[1:n],function(x){!isTRUE(all.equal(x,0))}))
	    vdiff2<-((s1[vtokeep]-s2[vtokeep])/((s1[vtokeep]+s2[vtokeep])))^2
	    if ((method=="logweighted_node_heights")&&(length(vtokeep)>2)){
		vdiff2[3:length(vtokeep)]<-vdiff2[3:length(vtokeep)]/log(3:length(vtokeep))
	    }
	    tree_dist<-sqrt(mean(vdiff2))
	}
	if (method=="bdcoeffs"){
	    if (requireNamespace("geiger",quietly=TRUE)){ 
    		tree2_avg_rate<-geiger::bd.km(tree2)
		tree1_avg_rate<-geiger::bd.km(tree1)
	    }else{
	    ## geiger is currently orphaned on CRAN so the geiger::bd.km() function was copied into pcmabc
    		tree2_avg_rate<-.geiger_bd.km(tree2)
		tree1_avg_rate<-.geiger_bd.km(tree1)	    
	    }
	## total variation distance between two exponentials -> average rate
	    L1<-max(c(tree1_avg_rate,tree2_avg_rate))
	    L2<-min(c(tree1_avg_rate,tree2_avg_rate))
	    if (L1!=L2){tree_dist<-0.5*((L2/L1)^(L2/(L1-L2)))*((L1-L2)/L1)}
	}

	if (method=="RF.dist"){
	    tree_dist<-phangorn::RF.dist(tree1,tree2)
	}
	if (method=="wRF.dist"){
	    tree_dist<-phangorn::wRF.dist(tree1,tree2,normalize=FALSE)
	}
	if (method=="wRFnorm.dist"){
	    tree_dist<-phangorn::wRF.dist(tree1,tree2,normalize=TRUE)
	}
	if (method=="KF.dist"){
	    tree1$tip.label<-tree2$tip.label[1:length(tree1$tip.label)]
	    tree_dist<-phangorn::KF.dist(tree1,tree2)
	}
	if (method=="path.dist"){
	    tree_dist<-phangorn::path.dist(tree1,tree2,use.weight=FALSE)
	}
	if (method=="path.dist.weights"){
	    tree_dist<-phangorn::path.dist(tree1,tree2,use.weight=TRUE)
	}
	if (method=="dist.topo.KF1994"){
	    tree_dist<-ape::dist.topo(tree1,tree2,method="score")
	}
    }
    tree_dist
}

.calc_data_dist<-function(data1,data2,method="order"){
## truedata is assumed to be data2, simulated data1
    data_dist<-0
    if (method=="order"){
	m<-ncol(data1)
	data1<-apply(data1,2,function(x){rev(x[order(abs(x))])})
	data1<-matrix(data1,ncol=m)

	m<-ncol(data2)
	data2<-apply(data2,2,function(x){rev(x[order(abs(x))])})
	data2<-matrix(data2,ncol=m)
		
	n<-min(nrow(data1),nrow(data2))
	data1<-c(apply(data1[1:n,,drop=FALSE],2,function(x){x/max(abs(x),na.rm=TRUE)}))
	data2<-c(apply(data2[1:n,,drop=FALSE],2,function(x){x/max(abs(x),na.rm=TRUE)}))	
	vtokeep<-which(sapply(abs(data1)+abs(data2),function(x){!isTRUE(all.equal(x,0))&&!is.na(x)&&!is.infinite(x)&&!is.nan(x)}))
	data1<-data1[vtokeep]
	data2<-data2[vtokeep]
	if ((length(data1)>0) && (length(data2)>0)){
	    vdiff2<-((data1-data2)/(abs(data1)+abs(data2)))^2
	    data_dist<-sqrt(mean(vdiff2,na.rm=TRUE))/2
	}else{
	    data_dist<-NA
	    if ((length(data1)==0) && (length(data2)==0)){data_dist<-0}
	}
    }
    if ((method=="variance")||(method=="variancemean")){
	if (ncol(data1)==1){mC1<-matrix(stats::var(data1[,1],na.rm=TRUE),1,1)}else{mC1<-stats::cov(data1,use="pairwise.complete.obs")}
	if (ncol(data2)==1){mC2<-matrix(stats::var(data2[,1],na.rm=TRUE),1,1)}else{mC2<-stats::cov(data2,use="pairwise.complete.obs")}

	p<-min(nrow(mC1),nrow(mC2))
	vC1<-c(mC1[1:p,1:p])
	vC2<-c(mC2[1:p,1:p])
	vtokeep<-which(sapply(abs(vC1)+abs(vC2),function(x){!isTRUE(all.equal(x,0))}))
	vC1<-vC1[vtokeep]
	vC2<-vC2[vtokeep]
	vdiff2<-((vC1-vC2)/(abs(vC1)+abs(vC2)))^2
	data_dist<-sqrt(mean(vdiff2))/2	
    }
    if (method=="variancemean"){
	vm1<-apply(data1,2,base::mean,na.rm=TRUE)
	vm2<-apply(data2,2,base::mean,na.rm=TRUE)
	p<-min(length(vm1),length(vm2))
	vm1<-c(vm1[1:p])
	vm2<-c(vm2[1:p])
	vtokeep<-which(sapply(abs(vm1)+abs(vm2),function(x){!isTRUE(all.equal(x,0))}))
	vm1<-vm1[vtokeep]
	vm2<-vm2[vtokeep]
	vdiff2<-((vm1-vm2)/(abs(vm1)+abs(vm2)))^2
	data_dist<-0.5*(data_dist+sqrt(mean(vdiff2))/2)
    }

    if (method=="bins"){
        truedata<-data2
        simdata<-data1
        numbins<-2*as.integer(nrow(truedata)/30)+5
	boundries<-apply(rbind(simdata,truedata),2,function(x,numbins){seq(from=min(x,na.rm=TRUE)-0.5,to=max(x,na.rm=TRUE)+0.5,length=numbins)},numbins=numbins)
	if (ncol(simdata)==1){boundries<-list(boundries)}
	if (!is.list(boundries)){
	    boundries<-sapply(1:ncol(boundries),function(i,data){data[,i]},data=boundries,simplify=FALSE)	
	}
	simhist<-sapply(1:ncol(simdata),function(i,x,boundries){graphics::hist(x[,i],plot=FALSE,breaks=boundries[[i]])$counts/length(x[,i])},x=simdata,boundries=boundries,simplify=FALSE)
	truehist<-sapply(1:ncol(truedata),function(i,x,boundries){graphics::hist(x[,i],plot=FALSE,breaks=boundries[[i]])$counts/length(x[,i])},x=truedata,boundries=boundries,simplify=FALSE)

	if (is.null(weights)){weights<-rep(1/length(simhist),length(simhist))}
	data_dist<-0.5*sum(sapply(1:length(simhist),function(i,x1,x2,w){w[i]*sum(abs(x1[[i]]-x2[[i]]))},x1=simhist,x2=truehist,w=weights,simplify=TRUE))

    }
    data_dist
}

.calc_ABCdist<-function(simtree,truetree,simdata,truedata,dist_method=c("variancemean","node_heights"),distfunc=.calc_separate_dist){
	#save(simtree,truetree,simdata,truedata,dist_method,file="qpa.RData");stop()
	distfunc(simtree,truetree,simdata,truedata,dist_method)
}

.calc_separate_dist<-function(simtree,truetree,simdata,truedata,dist_method=c("variancemean","node_heights")){
	dist<-(.calc_phyl_dist(simtree,truetree,method=dist_method[2])+.calc_data_dist(simdata,truedata,method=dist_method[1]))
	if (is.na(dist) || is.nan(dist) || is.infinite(dist)){dist<-100000} #normally the distance is scaled to lie in [0,1]
	else{if (!is.na(dist_method[2])){dist<-dist/2}}
	dist
}
