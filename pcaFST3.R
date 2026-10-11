TEST:
library(data.table)
a1f<-list.files(pattern="ad1_fff") 
a1f<-a1f[1:23] 
a2f<-a1f 
a2f<-gsub("ad1","ad2",a2f) 
N<-length(a1f) 
ids<-read.table("IDs.txt",header=FALSE) 
temp<-gsub("ad1_fff_o_lycpool_chrom","",a1f) 
chrom<-gsub(".txt","",temp) 
reps<-grep(pattern="rep",ids[,1])

## get rid of MEN and replicates 
reps <- grep("rep", ids[, 1])
men  <- grep("MEN", ids[, 1])
drop <- unique(c(reps, men))
keep <- setdiff(seq_len(nrow(ids)), drop)
print(ids[drop, 1])   # check: should list the replicates and Lyc-MEN12
pdf("WG_LG_PCAs.pdf",width=7,height=10.5) 
par(mfrow=c(3,2)) 
par(mar=c(4.5,5.5,2.5,1.5)) 
cl<-1.3;ca<-1.1;cm<-1.4 
pch_all <- rep(15:20, each = 7, length.out = nrow(ids))
pchs <- pch_all[keep]
cols <- ids[keep, 2]
labs <- ids[keep, 1]

## one PCA (PC1 vs PC2) per chromosome, 6 chromosomes per page
pdf("WG_LG_PCAs_noMEN.pdf", width = 7, height = 10.5)
par(mfrow = c(3, 2))
par(mar = c(4.5, 5.5, 2.5, 1.5))
cl <- 1.3; ca <- 1.1; cm <- 1.4

for(i in 1:N){ 
  a1<-as.matrix(fread(a1f[i],header=F)) 
  a2<-as.matrix(fread(a2f[i],header=F)) 
  n<-a1+a2 
  p<-a2/(a1+a2) ## non-ref
  p[is.na(p)]<-0.001  ## could this fuck with nusuem samples? Consider: p[is.na(p)] <- rowMeans(p, na.rm = TRUE)[row(p)[is.na(p)]]
  pc<-prcomp(t(p[,-c(reps,17)]),center=TRUE,scale=FALSE)
  pct <- round(summary(pc)$importance[2, 1:2] * 100, 1)	
	
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





Working ONE:


library(data.table)

a1f <- list.files(pattern = "ad1_fff")
a1f <- a1f[1:23]
a2f <- gsub("ad1", "ad2", a1f)
N <- length(a1f)
ids <- read.table("IDs.txt", header = FALSE)
temp <- gsub("ad1_fff_o_lycpool_chrom", "", a1f)
chrom <- gsub(".txt", "", temp)

## get rid of MEN and replicates
reps <- grep("rep", ids[, 1])
men  <- grep("MEN", ids[, 1])
drop <- unique(c(reps, men))
keep <- setdiff(seq_len(nrow(ids)), drop)
print(ids[drop, 1])   # check: should list the replicates and Lyc-MEN12

## symbols assigned on the full set, then subset, so each pop keeps its symbol
pch_all <- rep(15:20, each = 7, length.out = nrow(ids))
pchs <- pch_all[keep]
cols <- ids[keep, 2]
labs <- ids[keep, 1]

## one PCA (PC1 vs PC2) per chromosome, 6 chromosomes per page
pdf("WG_LG_PCAs_noMEN.pdf", width = 7, height = 10.5)
par(mfrow = c(3, 2))
par(mar = c(4.5, 5.5, 2.5, 1.5))
cl <- 1.3; ca <- 1.1; cm <- 1.4

for (i in 1:N) {
  a1 <- as.matrix(fread(a1f[i], header = FALSE))
  a2 <- as.matrix(fread(a2f[i], header = FALSE))
  p <- a2 / (a1 + a2)  ## non-ref allele frequency

  ## missing data: 0.001 treats "no data" as "fixed for ref" and can pull
  ## low-coverage (museum) pools together. Mean imputation is the safer option:
  # p[is.na(p)] <- rowMeans(p, na.rm = TRUE)[row(p)[is.na(p)]]
  p[is.na(p)] <- 0.001

  ## pca (replicates and MEN removed)
  pc <- prcomp(t(p[, keep]), center = TRUE, scale = FALSE)
  pct <- round(summary(pc)$importance[2, 1:2] * 100, 1)

  ## PC1 vs PC2
  plot(pc$x[, 1], pc$x[, 2], pch = pchs, col = cols,
       xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
       cex.lab = cl, cex.axis = ca)
  title(main = paste0("Chromosome ", chrom[i]), cex.main = cm)
  text(pc$x[, 1], pc$x[, 2], labs, cex = .7)
}
dev.off()
