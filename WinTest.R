## Zach's WinTest.R logic: average 5 consecutive windows, normalize the
## three topology scores to sum to 1, genome-wide means, stacked barplot.
## Run from inside a win<P1>x<P2>x<P3> folder:  Rscript ../WinPlot.R
args <- commandArgs(TRUE)
q    <- if (length(args)) args[1] else sub("^win", "", basename(getwd()))
tax  <- strsplit(q, "x")[[1]]
nwin <- 1
cols <- c("A+B", "A+C", "B+C")
labs <- c(paste0(tax[1], "+", tax[2]), paste0(tax[1], "+", tax[3]), paste0(tax[2], "+", tax[3]))

res <- list(); nneg <- 0
for (ch in 1:23) {
  d <- read.table(paste0("winout", ch, ".tsv"), header = TRUE, check.names = FALSE)
  s <- as.matrix(d[, cols])
  g <- ceiling(seq_len(nrow(s)) / nwin)
  avg <- rowsum(s, g) / as.vector(table(g))
  nneg <- nneg + sum(avg < 0)
  avg[avg < 0] <- 0
  p <- avg / rowSums(avg)
  res[[ch]] <- data.frame(chrom = ch, block = seq_len(nrow(p)),
                          start_col = as.vector(tapply(d$pos, g, min)),
                          p, check.names = FALSE)
}
all <- do.call(rbind, res)
all <- all[complete.cases(all), ]
cat("negative averaged scores set to 0:", nneg, "\n")
write.table(all, paste0("WinSummary_", q, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

gw <- colMeans(all[, cols]); names(gw) <- labs
cat("genome-wide mean proportions:\n"); print(round(gw, 3))

pdf(paste0("WinPlot_", q, ".pdf"), width = 14, height = 4)
newchr <- c(TRUE, diff(all$chrom) != 0)
bp <- barplot(t(as.matrix(all[, cols])), space = ifelse(newchr, 3, 0), border = NA, axisnames = FALSE,
              col = c("#1b9e77", "#d95f02", "#7570b3"), ylab = "Proportion",
              xlab = "Chromosome", axes = TRUE)
mids <- tapply(bp, all$chrom, mean)
axis(1, at = mids, labels = names(mids), tick = FALSE, cex.axis = 0.8)
legend("topright", inset = c(0, -0.18), xpd = TRUE, horiz = TRUE, bty = "n",
       fill = c("#1b9e77", "#d95f02", "#7570b3"), legend = labs, cex = 0.8)
dev.off()
cat("wrote WinPlot_", q, ".pdf and WinSummary_", q, ".tsv\n", sep = "")
