## This file is part of pcmabc

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at bartoszekkj@gmail.com .



.f_phenotype_time<-function(Tmax,X0,step){
    res<-rbind(seq(from=0,to=Tmax,by=step),seq(from=X0,to=X0+Tmax,by=step))
    if (is.vector(res)){res<-matrix(res,nrow=2,ncol=1)}
    res
}

draw_phylproc<-function(simulobj){
    tree.ape<-simulobj$tree
    phenotype.trajectory<-simulobj$phenotype
    root.branch.phenotype<-simulobj$root.branch.phenotype
    PhylTraitProcess<-vector("list",1)
    names(PhylTraitProcess)<-"FullTrajectory"
    tree.height<-tree.ape$tree.height
    if (length(phenotype.trajectory)>0){
	PhylTraitProcess$FullTrajectory<-sapply(1:length(phenotype.trajectory),function(i,phenotype.trajectory,tree.ape,tree.height){
					mData<-t(phenotype.trajectory[[i]])
					mData[,1]<-mData[,1]+(tree.height-tree.ape$node.heights[which(tree.ape$node.heights[,1]==tree.ape$edge[i,1]),2])				    
					list(trajectory=mData)
					},phenotype.trajectory=phenotype.trajectory,tree.ape=tree.ape,tree.height=tree.height,simplify=FALSE)
    }
    PhylTraitProcess$FullTrajectory[[length(phenotype.trajectory)+1]]<-list(trajectory=t(root.branch.phenotype))
    utils::capture.output(mvSLOUCH::drawPhylProcess(PhylTraitProcess,plotlayout=c(1,nrow(simulobj$phenotype[[1]])-1)))
    NA
}



