RCC<-function(type="nex",path=tk_choose.files(caption = "Choose the tree file"), output="resultsRCC.tre") {

# Wilkinson RCC

# librairies

	if(require(ape)) {library(ape)} else {install.packages("ape");library(ape)}
	if(require(phylobase)) {library(phylobase)} else {install.packages("phylobase");library(phylobase)}
	if(require(phytools)) {library(phytools)} else {install.packages("phytools");library(phytools)}
	if(require(tcltk)) {library(tcltk)} else {install.packages("tcltk");library(tcltk)}
	if(require(compare)) {library(compare)} else {install.packages("compare");library(compare)}

# Load trees

	path<-tk_choose.files(caption = "Choose the tree file")
	asP4<- function (x) {as(x,"phylo4")}
	if(type=="nex") {
		trees<-read.nexus(path)
	}
	if(type=="tre") {
		trees<-read.tree(path)
	}
	#for (i in 1:length(trees)){
	#	trees[[i]]$edge.length<-compute.brlen(trees[[i]], NA)
	#}
	trees<-lapply(trees,asP4)

# Components

	tree_index<-1
	hc_set <- "kiki"
	while (hc_set != 0) {
		incomp_fun  <- function (tree_cur,node_comp) {names(descendants(tree_cur, node_comp, type = "tips"))}
		outcomp_fun <- function (tree_cur,node_comp) {setdiff(tipLabels(tree_cur),names(descendants(tree_cur, node_comp, type = "tips")))}
		if (tree_index==1) {
			tree_1<-trees[[tree_index]]
			tree_2<-trees[[tree_index+1]]
			comp_tree_1<-setdiff(nodeId(tree_1,type="internal"),nodeId(tree_1,type="root"))
			comp_tree_2<-setdiff(nodeId(tree_2,type="internal"),nodeId(tree_2,type="root"))
			# incomp<-as.list(numeric(length(comp_tree_1)*length(comp_tree_2)));dim(incomp) <- c(length(comp_tree_1),length(comp_tree_2));rownames(incomp)<-comp_tree_1;colnames(incomp)<-comp_tree_2
			# find common taxa of inside sets ("in" hereafter) from trees 1 and 2
			combinTmp <- expand.grid(1:length(comp_tree_1), 1:length(comp_tree_2))
			incomp1<-list(); for (i in 1:length(comp_tree_1)) {incomp1[[i]]<-incomp_fun(tree_1,comp_tree_1[i])}; names(incomp1)<-comp_tree_1
			incomp2<-list(); for (i in 1:length(comp_tree_2)) {incomp2[[i]]<-incomp_fun(tree_2,comp_tree_2[i])}; names(incomp2)<-comp_tree_2
			outcomp1<-list(); for (i in 1:length(comp_tree_1)) {outcomp1[[i]]<-outcomp_fun(tree_1,comp_tree_1[i])}; names(outcomp1)<-comp_tree_1
			outcomp2<-list(); for (i in 1:length(comp_tree_2)) {outcomp2[[i]]<-outcomp_fun(tree_2,comp_tree_2[i])}; names(outcomp2)<-comp_tree_2
			taxin<-list()
			for (i in 1:nrow(combinTmp)) {
				taxin[[i]]<-intersect(incomp1[[combinTmp[i,1]]],incomp2[[combinTmp[i,2]]])
				if ((i%%1000)==0) {cat(paste("comb n° ",i," over ",nrow(combinTmp),"\n",sep=""))}
			}
			# find common taxa of outside sets ("out" hereafter) from trees 1 and 2
			taxout<-list()
			for (i in 1:nrow(combinTmp)) {
				taxout[[i]]<-intersect(outcomp1[[combinTmp[i,1]]],outcomp2[[combinTmp[i,2]]])
				if ((i%%1000)==0) {cat(paste("comb n° ",i," over ",nrow(combinTmp),"\n",sep=""))}
			}
			# conditions
			condition_in<-sapply(taxin, function (x) length(x)>=2)
			condition_out<-sapply(taxout, function (x) length(x)>=1)
			taxin<-taxin[condition_in&condition_out]
			taxout<-taxout[condition_in&condition_out]
			if(length(taxin) > 1) {cond<-TRUE} else {cond<-FALSE}
			while (cond) {				# while there is still pairs of identical nts
				nts_comb<-combn(c(1:length(taxin)),2)
				condition_red1<-logical(ncol(nts_comb)) # check that a nts x is equal to a nts y
				condition_red2<-logical(ncol(nts_comb)) # check that a nts y is equal to a nts x
				for(i in 1:ncol(nts_comb)) {
					condition_red1[i]<-length(intersect(unlist(taxin[nts_comb[1,i]]),unlist(taxin[nts_comb[2,i]]))) == length (unlist(taxin[nts_comb[1,i]])) & length(intersect(unlist(taxout[nts_comb[1,i]]),unlist(taxout[nts_comb[2,i]]))) == length (unlist(taxout[nts_comb[1,i]]))
					condition_red2[i]<-length(intersect(unlist(taxin[nts_comb[1,i]]),unlist(taxin[nts_comb[2,i]]))) == length (unlist(taxin[nts_comb[2,i]])) & length(intersect(unlist(taxout[nts_comb[1,i]]),unlist(taxout[nts_comb[2,i]]))) == length (unlist(taxout[nts_comb[2,i]]))
					if ( condition_red1[i]==T && condition_red2[i]==T ) condition_red1[i]<-FALSE
					if ((i%%100000)==0) {cat(paste("comb ",i," over ",ncol(nts_comb), " pairs of nts\n",sep=""))}
				}
				redund1<-unique(nts_comb[1,which(condition_red1)])
				redund2<-unique(nts_comb[2,which(condition_red2)])
				taxin<-taxin[setdiff(c(1:length(taxin)), sort(unique(c(redund1,redund2))))]
				taxout<-taxout[setdiff(c(1:length(taxout)), sort(unique(c(redund1,redund2))))]				
				if(all(!c(condition_red1,condition_red2)) | length(taxin) == 1) cond<-FALSE
			}
			tree_index<-tree_index+1 ; cat(paste("tree index = ",tree_index,"\n",sep=""))
			if(tree_index == length(trees)) hc_set<-0
		} else {
			tree_n<-trees[[tree_index+1]]
			comp_tree_n<-setdiff(nodeId(tree_n,type="internal"),nodeId(tree_n,type="root"))
			combinTmp <- expand.grid(1:length(taxin), 1:length(comp_tree_n))
			incompn<-list(); for (i in 1:length(comp_tree_n)) {incompn[[i]]<-incomp_fun(tree_n,comp_tree_n[i])}; names(incompn)<-comp_tree_n
			outcompn<-list(); for (i in 1:length(comp_tree_n)) {outcompn[[i]]<-outcomp_fun(tree_n,comp_tree_n[i])}; names(outcompn)<-comp_tree_n
			taxintemp<-list()
			for (i in 1:nrow(combinTmp)) {
				taxintemp[[i]]<-intersect(incompn[[combinTmp[i,2]]],taxin[[combinTmp[i,1]]])
				if ((i%%1000)==0) {cat(paste("comb n° ",i," over ",nrow(combinTmp),"\n",sep=""))}
			}
			# trouver les taxons commun ÃƒÂ  l'ext des composantes des arbres 1 et 2
			taxouttemp<-list()
			for (i in 1:nrow(combinTmp)) {
				taxouttemp[[i]]<-intersect(outcompn[[combinTmp[i,2]]],taxout[[combinTmp[i,1]]])
				if ((i%%1000)==0) {cat(paste("comb n° ",i," over ",nrow(combinTmp),"\n",sep=""))}
			}
			taxout<-taxouttemp
			taxin<-taxintemp
			# conditions
			condition_in<-sapply(taxin, function (x) length(x)>=2)
			condition_out<-sapply(taxout, function (x) length(x)>=1)
			taxin<-taxin[condition_in&condition_out]
			taxout<-taxout[condition_in&condition_out]
			try(if(length(taxin)==0) stop("\n\nUnfortunately, even the RCC cannot find any signal. You have a star phylogeny.\n\n"))
			if(length(taxin) > 1) {cond<-TRUE} else {cond<-FALSE}
			while (cond) {				# while there is still pairs of identical nts
				nts_comb<-combn(c(1:length(taxin)),2)
				condition_red1<-logical(ncol(nts_comb)) # check that a nts x is equal to a nts y
				condition_red2<-logical(ncol(nts_comb)) # check that a nts y is equal to a nts x
				for(i in 1:ncol(nts_comb)) {
					condition_red1[i]<-length(intersect(unlist(taxin[nts_comb[1,i]]),unlist(taxin[nts_comb[2,i]]))) == length (unlist(taxin[nts_comb[1,i]])) & length(intersect(unlist(taxout[nts_comb[1,i]]),unlist(taxout[nts_comb[2,i]]))) == length (unlist(taxout[nts_comb[1,i]]))
					condition_red2[i]<-length(intersect(unlist(taxin[nts_comb[1,i]]),unlist(taxin[nts_comb[2,i]]))) == length (unlist(taxin[nts_comb[2,i]])) & length(intersect(unlist(taxout[nts_comb[1,i]]),unlist(taxout[nts_comb[2,i]]))) == length (unlist(taxout[nts_comb[2,i]]))
					if ( condition_red1[i]==T && condition_red2[i]==T ) condition_red1[i]<-FALSE
					if ((i%%100000)==0) {cat(paste("comb ",i," over ",ncol(nts_comb), " pairs of nts\n",sep=""))}
				}
				redund1<-unique(nts_comb[1,which(condition_red1)])
				redund2<-unique(nts_comb[2,which(condition_red2)])
				taxin<-taxin[setdiff(c(1:length(taxin)), sort(unique(c(redund1,redund2))))]
				taxout<-taxout[setdiff(c(1:length(taxout)), sort(unique(c(redund1,redund2))))]				
				if(all(!c(condition_red1,condition_red2)) | length(taxin) == 1) cond<-FALSE
			}
			if(tree_index == length(trees)-1) hc_set<-0
			tree_index<-tree_index+1; cat(paste("tree index = ",tree_index,"\n",sep=""))
		}
	}


# set compatibility
	taxinF<-taxin;taxoutF<-taxout
	RCCprofile<-list()
	lists<-list()
	lists<-Map(c, taxoutF, taxinF)
	cond<-TRUE
	while (cond & length(lists)>1) {				# while there is still pairs of identical lists
		nts_comb<-combn(c(1:length(lists)),2)
		condition_red1<-logical(ncol(nts_comb))
		for(i in 1:ncol(nts_comb)) {
			condition_red1[i]<-compare(lists[nts_comb[1,i]], lists[nts_comb[2,i]], allowAll = TRUE, shorten = FALSE)$result
		}
		nts_comb[,condition_red1]
		redund<-unique(nts_comb[1,which(condition_red1)])
		lists<-lists[setdiff(c(1:length(lists)), unique(redund))]
		if(all(!condition_red1)) cond<-FALSE
	}
	#lists <- apply( cbind(taxoutF, taxinF) , 1 , unlist) # complete sets
	lists <- unique(lapply(lists,unname))	# unique complete sets (remove duplicated ones)
	for (j in 1:length(lists)) {
		RCCtempin<-RCCtempout<-list()
		RCCprofile[[j]]<-starTree(lists[[j]])
		for (i in 1:length(taxinF)) {
			if( length(intersect(c(taxoutF[[i]],taxinF[[i]]),lists[[j]])) == length(lists[[j]]) ) {
				RCCtempin[[i]]<-intersect(taxinF[[i]],lists[[j]])
				RCCtempout[[i]]<-intersect(taxoutF[[i]],lists[[j]])
				if(length(RCCtempin[[i]])>=2 & length(RCCtempout[[i]])>=1) {
					TEdge<-RCCprofile[[j]]$edge
					TLabel<-RCCprofile[[j]]$tip.label
					newanc<-unique(TEdge[match(RCCtempin[[i]],TLabel),1])
					TEdge[match(RCCtempin[[i]],TLabel),1]<-max(TEdge[,1])+1
					TEdge<-rbind(TEdge,c(newanc,max(TEdge[,1])))
					RCCprofile[[j]]$edge<-TEdge
					RCCprofile[[j]]$Nnode<-RCCprofile[[j]]$Nnode+1
				}
			}
		}
	}

# plot Consensus + RCC profile
	# absent taxa
	absent<-list()
	for (i in 1:length(RCCprofile)) {
		absent[[i]]<-setdiff(trees[[1]]@label, RCCprofile[[i]]$tip.label)
	}
	Absentlengths<-NA;for(i in 1:length(absent)) {Absentlengths[i]<-length(absent[[i]])}
	Absentorder<-cbind(Absentlengths,c(1:length(absent)))
	Absentindex<-Absentorder[order(Absentorder[,1],decreasing=FALSE),]
	if(!is.null(dim(Absentindex))) absent<-absent[Absentindex[,2]]
	# strict consensus
	asP<- function (x) {as(x,"phylo")}
	cons<-consensus(lapply(trees,asP))
	# Order by trees lengths 
	RCClengths<-NA;for(i in 1:length(RCCprofile)) {RCClengths[i]<-length(RCCprofile[[i]]$tip.label)}
	RCCorder<-cbind(RCClengths,c(1:length(RCCprofile)))
	RCCindex<-RCCorder[order(RCCorder[,1],decreasing=TRUE),]
	if(!is.null(dim(RCCindex))) RCCprofile<-RCCprofile[RCCindex[,2]]
	# add strict consensus if not in RCC profile
	if (cons$Nnode == 1) { 
		RCCprofile<-c(list(cons),RCCprofile)
		absent<-c(list(NA),absent)
	}
	#graphical parameters
	ntrees<-length(RCCprofile)
	pgdiv<-0
	if (max(RCClengths)>20) {pgdiv<-1}
	# plot trees
	t<-1
	if (pgdiv==1) {
		for(i in 1:length(RCCprofile)) {
			dev.new(width=20, height=16)
			par(mfrow=c(1,2))
			#par(fig=c(0,0.6,0,1),mar=c(0,0,1,0))
			plot(ladderize(RCCprofile[[t]]),cex=.5)
			if(t==1) {title("Strict consensus")} else {title(paste("RCC profile - tree",t-1))}
			#par(fig=c(0.6,1,0,1),mar=c(0,0,0,0), new=TRUE)
			plot(1, type="n", axes=F, xlab="", ylab="")
			if (length(absent[[t]])>0) {for (a in 1:length(absent[[t]])) {text(grconvertX(0.5,"npc"),grconvertX(a/(length(absent[[t]])+2),"npc"),absent[[t]][a],cex=.5)};text(grconvertX(0.5,"npc"),grconvertX(0.9,"npc"),"Absent:")}
			t<-t+1
			if(t!=length(ntrees)) x11()
		}
	} else {
		for(i in 1:ceiling(ntrees/2)) {
			for(j in 1:2){
				if(t<=length(RCCprofile)) {
					if (j==1) par(mfrow=c(2,2))
					#par(fig=c(0,0.6,0,1),mar=c(0,0,1,0))
					plot(ladderize(RCCprofile[[t]]))
					if(t==1) {title("Strict consensus")} else {title(paste("RCC profile - tree",t-1))}
					#par(fig=c(0.6,1,0,1),mar=c(0,0,0,0), new=TRUE)
					plot(1, type="n", axes=F, xlab="", ylab="")
					if (length(absent[[t]])>0) {
						for (a in 1:length(absent[[t]])) {
							text(grconvertX(0.5,"npc"),grconvertX(a/(length(absent[[t]])+2),"npc"),absent[[t]][a])
						}
						text(grconvertX(0.5,"npc"),grconvertX(0.9,"npc"),"Absent:")
					}
					if((t%%2)==0 & t!=ntrees)  x11()
					t<-t+1
				}
			}
		}
	}
	class(RCCprofile) <- "multiPhylo"
	if(length(RCCprofile)==1){names(RCCprofile)<-"Strict consensus "}else {names(RCCprofile) <- c("Strict consensus ", paste("RCC_",1:(length(RCCprofile)-1)," ",sep=""))}
	write.tree(RCCprofile,output, tree.names=TRUE)
	return(list(RCCprofile,absent))
}
