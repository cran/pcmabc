## This file is part of pcmabc

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at bartoszekkj@gmail.com .



get_phylogenetic_sample<-function(pcmabc_simulobj,bOnlyContemporary=FALSE,tol=.Machine$double.eps^0.5){
    
    vTipNodes<-NA
    if (bOnlyContemporary){
        vTipNodes<-which(sapply(pcmabc_simulobj[[1]]$node.heights,function(x,tol){isTRUE(all.equal(x,0,tolerance=tol))},tol=tol,simplify=TRUE))
    }else{
	vTipNodes<-setdiff(unique(pcmabc_simulobj[[1]]$edge[,2]),unique(pcmabc_simulobj[[1]]$edge[,1]))
    }
    vTipedges<-sapply(vTipNodes,function(x,mTreeEdge){which(mTreeEdge[,2]==x)},mTreeEdge=pcmabc_simulobj[[1]]$edge,simplify=TRUE)
    
    mContemporarySample<-t(sapply(vTipedges,function(i,lphenotype){
	lphenotype[[i]][-1,ncol(lphenotype[[i]])]
    },lphenotype=pcmabc_simulobj[[2]],simplify=TRUE))
    if (nrow(mContemporarySample)==1){mContemporarySample<-t(mContemporarySample)}
    mContemporarySample
}

