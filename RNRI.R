## Rescaled Nodal Retention Index

RNRI<-function(IZ=file.choose(),IA=file.choose(),size.dots=5) {

# crgelent librairies
if(require(ape)) {library(ape)} else {install.packages("ape",dependencies = TRUE);library(ape)}
if(require(phylobase)) {library(phylobase)} else {install.packages("phylobase",dependencies = TRUE);library(phylobase)}
if(require(phytools)) {library(phytools)} else {install.packages("phytools",dependencies = TRUE);library(phytools)}
if(require(tcltk)) {library(tcltk)} else {install.packages("tcltk",dependencies = TRUE);library(tcltk)}

# importer fichiers
	# import characters (3ia)
	IA<-readLines(tk_choose.files(caption = "Choose the characters file (.ia)"))
	tax <- 0 # booléen tant qu'il y a encore des taxons
	clineTax <- grep("Taxa",IA)
	numtaxlisb<-NA
	namestax<-NA
	ind<-1
	while (tax == 0) {
		if (!substr(IA[clineTax],1,1)==";") {
			clineTax<-clineTax+1
			numtaxlisb[ind]<-as.numeric(substr(IA[clineTax],3,regexpr(" = ",IA[clineTax])-1))
			namestax[ind]<-substr(IA[clineTax],regexpr(" = ",IA[clineTax])+3,nchar(IA[clineTax]))
			ind<-ind+1
		} else {tax<-1;numtaxlisb<-numtaxlisb[-length(numtaxlisb)];namestax<-namestax[-length(namestax)]}
	}
	corres<-cbind(numtaxlisb,namestax)
	cline <- grep("Characters",IA)
	ch <- 0 # booléen tant qu'il y a encore des compatibles
	chars<-list() # liste des arbres compatibles
	j<-1 # indice liste
	cline<-cline+1
	while (ch == 0) {
		char_temp <- IA[cline]
		if(min(regexpr("[(]",char_temp ))==-1){
			if(min(regexpr("[;]",char_temp ))!=-1){ ch<-1 } else {
				cline<-cline+1
			}
		} else {
			start<- min(regexpr("[(]",char_temp ))
			end<- max(gregexpr("[)]",char_temp )[[1]])
			chars[[j]]<-as(read.tree(text=paste(substr(char_temp,start,end),";",sep="")),"phylo4")
			labels(chars[[j]],"tip")<-corres[,2][as.numeric(labels(chars[[j]],"tip"))+1]
			j<-j+1
			cline<-cline+1
		}
	}

	# import arbre (3iz)
	path<-tk_choose.files(caption = "Choose the tree file")
	IZ<-readLines(path)
	if(length(grep("-<D06>-",IZ))>0) {
		taxA <- 0 # booléen tant qu'il y a encore des taxons
		clineTax <- grep("Taxa",IZ)+1
		numtaxlisb<-NA
		namestax<-NA
		ind<-1
		while (taxA == 0) {
			if (substr(IZ[clineTax],1,1)==".") {
				numtaxlisb[ind]<-as.numeric(substr(IZ[clineTax],20,21))
				namestax[ind]<-substr(IZ[clineTax],24,nchar(IZ[clineTax]))
				clineTax<-clineTax+1
				ind<-ind+1
			} else {taxA<-1}
		}
		corres<-cbind(numtaxlisb,namestax)
		cline <- grep("Compatible trees :",IZ)
		ct <- 0 # booléen tant qu'il y a encore des compatibles
		trees<-list() # liste des arbres compatibles
		i<-1 # indice liste
		while (ct == 0) {
			cline<-cline+1
			tree_temp <- IZ[cline]
			if (substr(tree_temp,1,1)==":") {
				start<- min(regexpr("[(]",tree_temp))
				end<- max(gregexpr("[)]",tree_temp)[[1]])
				trees[[i]]<-as(read.tree(text=paste(substr(tree_temp,start,end),";",sep="")),"phylo4")
				labels(trees[[i]],"tip")<-corres[,2][as.numeric(labels(trees[[i]],"tip"))+1]
				i<-i+1
			} else {ct<-1}
		}
	} else {
		asP4<- function (x) {as(x,"phylo4")}
		trees<-read.nexus(path)
		if(class(trees)=="phylo") {
			tree<-list()
			trees<-as(trees,"phylo4")
			tree[[1]]<-trees
			trees<-tree
		} else {
			trees<-lapply(trees,asP4)
		}
	}

# pour chaque arbre de compatible
tree_index<-1
restrees<-list()
for (i in 1:length(trees)) {
	tree_c<-trees[[i]]
	index_c<-1
	INDICES<-matrix(NA,length(setdiff(nodeId(tree_c,type="internal"),nodeId(tree_c,type="root"))),3)
	for (c in setdiff(nodeId(tree_c,type="internal"),nodeId(tree_c,type="root")) ){  # composante c
		INComp<-names(descendants(tree_c, c, type = "tips"))
		OUTComp<-setdiff(tipLabels(tree_c),INComp)
		INDICES[index_c,3]<-(length(INComp)-1)*(length(OUTComp))
		char_index<-char_numb<-1
		reschar<-matrix(NA,1,6);colnames(reschar)<-c("char index","char numb","tj","tcUj","pj","pcUj")
		for (char in chars) {			# chaque caractère
			#rootNode(char)<-nodeId(char,"internal")[is.na(ancestor(char,nodeId(char,"internal")))]
			char_c<-nodeId(char,"internal")[!is.na(ancestor(char,nodeId(char,"internal")))] # selection internal nodes
			for(s in char_c) {
				reschartemp<-rep(NA,6)
				INChar_c_s<-names(descendants(char, s, type = "tips"))
				OUTChar_c_s<-setdiff(tipLabels(char),INChar_c_s)
				reschartemp[1]<-char_index
				reschartemp[2]<-char_numb
				reschartemp[3]<-length(tipLabels(char))
				reschartemp[5]<-length(INChar_c_s)
				pcUj<-length(intersect(INComp,INChar_c_s))
				tcUj<-length(intersect(OUTComp,OUTChar_c_s))+pcUj
				reschartemp[4]<-tcUj
				if (pcUj==0) {reschartemp[6]<-1} else {reschartemp[6]<-pcUj}
				reschar<-rbind(reschar,reschartemp)
				char_index<-char_index+1
			}
			char_numb<-char_numb+1
		}
		reschar<-reschar[-1,]
		INDICES[index_c,1]<-c
		if(is.vector(reschar)) {INDICES[index_c,2]<-sum((reschar[6]-1)*(reschar[4]-reschar[6]))/sum((reschar[5]-1)*(reschar[3]-reschar[5]))
			} else {INDICES[index_c,2]<-sum((reschar[,6]-1)*(reschar[,4]-reschar[,6]))/sum((reschar[,5]-1)*(reschar[,3]-reschar[,5])) 
		}
		#edgeLength(trees[[i]],node=c)<-NRI
		index_c<-index_c+1
	}
	INFctotal<-sum(INDICES[,3])
	INDICES[,3]<-INDICES[,3]/INFctotal
	INDICES<-cbind(INDICES,INDICES[,2]/INDICES[,3])
	colnames(INDICES)<-c("node","NRIc","INFc","RNRIc")
	#INDICESs<-data.frame(nri=round(INDICES[,2],2),row.names=INDICES[,1]);INDICESs<-rbind(INDICESs,NA);rownames(INDICESs)[length(rownames(INDICESs))]<-setdiff(nodeId(tree_c,type="internal"),as.numeric(rownames(INDICESs)))
	#trees[[i]]<-addData(trees[[i]],node.data=INDICESs)
	#treePlot(trees[[i]])
	#
	#nodelabels(pie=range01(trees[[i]]@data),piecol=c("black","white"), cex=0.6)
	#edgeLength(trees[[i]])<-1
	#coor<-phyloXXYY(trees[[i]])
	tre<-as(trees[[i]],"phylo")
	x11();plot(ladderize(tre));title(i)
	#range01 <- function(x){(x-(min(x)))/((max(x)-min(x)))}
	range01 <- function(x){(x)/((max(x)-min(x)))}
	nodelabels(node=INDICES[,1],pie=rep(1,length(INDICES[,1])),cex=size.dots*range01(INDICES[,4]))
	nodelabels(round(INDICES[,4],2),node=INDICES[,1],adj=c(1,-0.2),frame="none")
	#text(coor$xx[1:length(nodeId(tree_c,type="internal"))]+0.05,coor$yy[1:length(nodeId(tree_c,type="internal"))],round(INDICES[,4],2))
	#plot(coor$xx+0.05,coor$yy,pch=20)
	tree_index<-tree_index+1
}
return(INDICES)
}
