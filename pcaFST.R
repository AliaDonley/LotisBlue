getwd() setwd()

library(data.table)
a1f<-list.files(pattern="ad1_fff") 
a1f<-a1f[1:23] 
a2f<-a1f 
a2f<-gsub("ad1","ad2",a2f) 
N<-length(a1f) 
ids<-read.table("IDs.txt",header=FALSE) 
temp<-gsub("ad1_fff_lycSpecPool_chrom","",a1f) 
chrom<-gsub(".txt","",temp) 
reps<-grep(pattern="rep",ids[,1])

pdf("WG_LG_PCAs.pdf",width=7,height=10.5) 
par(mfrow=c(3,2)) par(mar=c(4.5,5.5,2.5,1.5)) 
cl<-1.3;ca<-1.1;cm<-1.4 
for(i in 1:N){ 
  a1<-as.matrix(fread(a1f[i],header=F)) 
  a2<-as.matrix(fread(a2f[i],header=F)) 
  n<-a1+a2 
  p<-a2/(a1+a2) ## non-ref
  p[is.na(p)]<-0.001  ## could this fuck with nusuem samples? Consider: p[is.na(p)] <- rowMeans(p, na.rm = TRUE)[row(p)[is.na(p)]]
	## pca
	pc<-prcomp(t(p[,-c(reps,17)]),center=TRUE,scale=FALSE) ## drop replicates and MEN12
	o<-summary(pc)
	pct<-round(o$importance[2,1:3] * 100,1)

	plot(pc$x[,1],pc$x[,2],pch=rep(15:20,each=7),col=ids[-c(reps,17),2],xlab=paste("PC1 (",pct[1],")"),ylab=paste("PC2 (",pct[2],")"),cex.lab=cl,cex.axis=ca)
	title(main=paste("Chromosome ",chrom[i],sep=""),cex.main=cm)
	text(pc$x[,1],pc$x[,2],ids[-c(reps,17),1],cex=.7)
	xa<-min(pc$x[,1]) * .9
	ya<-pc$x[22,2] * .9
	#if(ya > 0){
	#	legend(xa,ya,ids[,1],pch=rep(15:20,each=7),col=ids[,2],ncol=3,cex=.6)
	#}
	plot(pc$x[,1],pc$x[,3],pch=rep(15:20,each=7),col=ids[-c(reps,17),2],xlab=paste("PC1 (",pct[1],")"),ylab=paste("PC3 (",pct[3],")"),cex.lab=cl,cex.axis=ca)
	title(main=paste("Chromosome ",chrom[i],sep=""),cex.main=cm)
	text(pc$x[,1],pc$x[,3],ids[-c(reps,17),1],cex=.7)

}
dev.off()



library(data.table)

a1f <- list.files(pattern = "ad1_fff")
a1f <- a1f[1:23]
a2f <- gsub("ad1", "ad2", a1f)
N <- length(a1f)
ids <- read.table("IDs.txt", header = FALSE)
temp <- gsub("ad1_fff_o_lycpool_chrom", "", a1f)
chrom <- gsub(".txt", "", temp)

## populations to drop: replicates + the MEN outgroup (matched by name, e.g. "Lyc-MEN12")
reps <- grep("rep", ids[, 1])
men  <- grep("MEN", ids[, 1])
drop <- unique(c(reps, men))
keep <- setdiff(seq_len(nrow(ids)), drop)
print(ids[drop, 1])   # check: should list the replicates and Lyc-MEN12

## assign symbols on the full set, then subset, so each pop keeps its symbol
pch_all <- rep(15:20, each = 7, length.out = nrow(ids))
pchs <- pch_all[keep]
cols <- ids[keep, 2]
labs <- ids[keep, 1]

pdf("WG_LG_PCAs2.pdf", width = 7, height = 10.5)
par(mfrow = c(3, 2))
par(mar = c(4.5, 5.5, 2.5, 1.5))
cl <- 1.3; ca <- 1.1; cm <- 1.4

for(i in 1:N){
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n<-a1+a2
	p<-a2/(a1+a2) ## non-ref

	p[is.na(p)]<-0.001
	## pca
	pc<-prcomp(t(p[,-c(reps,17)]),center=TRUE,scale=FALSE) ## drop replicates and MEN12
	o<-summary(pc)
	pct<-round(o$importance[2,1:3] * 100,1)

	plot(pc$x[,1],pc$x[,2],pch=rep(15:20,each=7),col=ids[-c(reps,17),2],xlab=paste("PC1 (",pct[1],")"),ylab=paste("PC2 (",pct[2],")"),cex.lab=cl,cex.axis=ca)
	title(main=paste("Chromosome ",chrom[i],sep=""),cex.main=cm)
	text(pc$x[,1],pc$x[,2],ids[-c(reps,17),1],cex=.7)
	xa<-min(pc$x[,1]) * .9
	ya<-pc$x[22,2] * .9
	#if(ya > 0){
	#	legend(xa,ya,ids[,1],pch=rep(15:20,each=7),col=ids[,2],ncol=3,cex=.6)
	#}
	plot(pc$x[,1],pc$x[,3],pch=rep(15:20,each=7),col=ids[-c(reps,17),2],xlab=paste("PC1 (",pct[1],")"),ylab=paste("PC3 (",pct[3],")"),cex.lab=cl,cex.axis=ca)
	title(main=paste("Chromosome ",chrom[i],sep=""),cex.main=cm)
	text(pc$x[,1],pc$x[,3],ids[-c(reps,17),1],cex=.7)

}
dev.off()
