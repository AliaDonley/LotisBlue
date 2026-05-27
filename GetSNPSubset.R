library(data.table)

miss <- vector("list", 23)
for(i in 1:23){
    ifile <- paste0("text_chrom_noBAT49_filteredV2_", i, ".fasta")
    dat <- fread(ifile, header=FALSE)
    miss[[i]] <- apply(dat==5, 2, mean)
}

# Total SNPs from filteredV2
total_snps <- sum(c(22277, 20489, 18858, 19902, 23617, 18437, 18266, 19136,
                    15003, 17950, 17563, 16639, 19164, 14663, 17616, 16497,
                    17742, 15549, 19042, 14712, 13759, 10907, 15198))
prop <- 6000 / total_snps
cat("prop =", prop, "\n")
cat("total_snps =", total_snps, "\n")

keepSNPs <- vector("list", 23)
for(i in 1:23){
    xx <- which(miss[[i]] == 0)
    n_keep <- floor(length(miss[[i]]) * prop)
    cat("Chrom", i, "- keeping:", n_keep, "from", length(xx), "complete SNPs\n")
    keepSNPs[[i]] <- sort(sample(xx, n_keep, replace=FALSE))
}

for(i in 1:23){
    out <- paste0("keepSNPs_maxProp_noBAT49_filteredV2_chrom", i)
    write.table(keepSNPs[[i]], file=out, row.names=FALSE, col.names=FALSE, quote=FALSE)
}

save(list=ls(), file="snps_maxProp_noBAT49_filteredV2.rdat")
