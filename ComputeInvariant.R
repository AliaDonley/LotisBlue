# Step 1: read genome base counts
cnts <- read.table("baseCounts_noBAT49_V2.txt", header=FALSE)
dim(cnts)

# Step 2: get 23 big chromosomes and sum base counts
totals <- apply(cnts[,-1], 1, sum)
chr <- which(totals >= 9211676)
cat("Number of chromosomes:", length(chr), "\n")
bcnt <- apply(cnts[chr,-1], 2, sum)
cat("Total genome base counts (A, C, G, T):", bcnt, "\n")

# Step 3: scale by proportion
prop <- 5989 / 402986  # noBAT49_noTBY51 total SNPs
cat("Proportion:", prop, "\n")
sbcnt <- floor(bcnt * prop)
cat("Scaled genome counts (A, C, G, T):", sbcnt, "\n")

# Step 4: read SNP counts
snps <- read.table("snpCounts_noBAT49_noTBY51.txt", header=FALSE)
dim(snps)  # should be 25 x 5
snpCnts <- floor(apply(snps[,-1], 2, mean))
cat("Mean SNP counts (A, C, G, T):", snpCnts, "\n")

# Step 5: calculate invariant counts
invar <- sbcnt - snpCnts
cat("Invariant counts (A, C, G, T):", invar, "\n")
