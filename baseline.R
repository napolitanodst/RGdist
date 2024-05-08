#Baseline

counts_matrix <- assays(spe)$counts
zero_x_gene <- c()
for(i in 1:nrow(counts_matrix)){
  zero_x_gene[i] = sum(counts_matrix[i,] == 0)
}

no_zero_x_gene <- 2310 - zero_x_gene
rowData(spe)$no_zero_x_gene <- no_zero_x_gene

# before
logcounts_matrix <- assays(spe)$logcounts
logcounts_matrix <- as.matrix(logcounts_matrix)
baseline_before <- rowSums(logcounts_matrix)/no_zero_x_gene
mean_expr_before <- mean(baseline_before)

dist_mean_before <- c()
for(i in 1:length(baseline_before)){
  dist_mean_before[i] = abs(baseline_before[i] - mean_expr_before)
}

rank_baseline_before <- rank(dist_mean_before, ties.method = "min")
rank_baseline_before <- 16187 - rank_baseline_before

rowData(spe)$baseline_before <- baseline_before
rowData(spe)$rank_baseline_before <- rank_baseline_before

# after
BClogcounts <- assays(spe)$BClogcounts
BClogcounts <- as.matrix(BClogcounts)
baseline_after <- rowSums(BClogcounts)/no_zero_x_gene
mean_expr_after <- mean(baseline_after)

dist_mean_after <- c()
for(i in 1:length(baseline_after)){
  dist_mean_after[i] = abs(baseline_after[i] - mean_expr_after)
}

rank_baseline_after <- rank(dist_mean_after, ties.method = "min")
rank_baseline_after <- 16187 - rank_baseline_after

rowData(spe)$baseline_after <- baseline_after
rowData(spe)$rank_baseline_after <- rank_baseline_after

baseline <- as.matrix(rowData(spe)[,c(3,4)])

baseline$baseline_before <- baseline_before
baseline$rank_baseline_before <- rank_baseline_before
baseline$baseline_after <- baseline_after
baseline$rank_baseline_after <- rank_baseline_after
