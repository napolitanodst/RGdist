# SVGs Identification With SpatialDE
# https://www.bioconductor.org/packages/release/bioc/vignettes/spatialDE/inst/doc/spatialDE.html

library(SpatialExperiment)
library(spatialDE)
library(ggspavis)

## BEFORE BATCH CORRECTION

#Prepare input data
sample_info <- as.data.frame(matrix(data = NA, 
                                    nrow = ncol(spe), 
                                    ncol = 3,
                                    dimnames = list(c(colnames(spe)),
                                                    c("x", "y", "total_counts"))))
sample_info[,"x"] <- spe$array_row
sample_info[,"y"] <- spe$array_col
expr_counts <- as.matrix(assay(spe, "counts"))
sample_info[, "total_counts"] <- colSums(expr_counts)
coordinates <- sample_info[, c("x", "y")]

#Stabilize
expr_norm <- stabilize(expr_counts)

#Regress Out
expr_resid <- regress_out(expr_norm, sample_info = sample_info)

#Run
set.seed(123)
results <- spatialDE::run(expr_resid, coordinates = coordinates)


## AFTER BATCH CORRECTION

#Prepare input data
expr_counts_bc <- as.matrix(assay(spe, "BCcounts"))
sample_info[, "total_counts"] <- colSums(expr_counts_bc)

#Stabilize
expr_norm_bc <- stabilize(expr_counts_bc)

#Regress Out
expr_resid_bc <- regress_out(expr_norm_bc, sample_info = sample_info)

#Run
set.seed(123)
results_bc <- spatialDE::run(expr_resid_bc, coordinates = coordinates)

# Ranking
results$rank <- rank(results$qval, ties.method = "min")
results_bc$rank <- rank(results2$qval, ties.method = "min")

save(results, results_bc, file = "spatialDE_results.RData")
