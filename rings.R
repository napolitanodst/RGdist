
## requirements:
# BiocManager::install("SpatialExperiment")

config <- read.table("config.txt")[,2]
names(config) <- read.table("config.txt")[,1]


library(SpatialExperiment)
load(file.path(config["input_path"], "spe-spots.RData"))

counts_matrix <- read.csv(file.path(config["input_path"], "counts_matrix.csv"), header = T, row.names = 1)

# How many rings?
hist(spots$distance,
     main = "Distances from injection site",
     xlab = "Distances",
     col = "darkolivegreen2")

#Looking at the distances histogram we decided to start using a bin = 200 through which we got 15 rings.

# 15 rings
rings <- .bincode(spots$distance, c(seq(from = -1, to = 2967, by = 200), Inf))
spots$rings <- rings


#First of all, ripetere tutte le analisi utilizzando come prova i 15 rings. Se i risultati sono buoni, bisognerà capire qual è l'ampiezza ottimale per i rings (non possiamo deciderlo arbitrariamente senza validare).

# La matrice è troppo grande per cui è necessario lavorare in parallelo
library(parallel)

# Con 10 cpu il calcolo si riduce a 1.62014 mins.
num_cores <- as.numeric(config["num_cores"])

# Divido gli indici dei geni in gruppi per parallelizzazione
gene_groups <- split(1:nrow(counts_matrix), rep(1:num_cores, each = ceiling(nrow(counts_matrix) / num_cores))[1:nrow(counts_matrix)])

gene_ring_expr_parallel <- function(counts_matrix_row, rings) {
  urings <- unique(rings)
  ring_expr <- rep(NA, length(urings))
  
  for (i in 1:length(urings)) {
    index_bin <- urings[i] == rings
    ring <- counts_matrix_row[index_bin]
    no_zero <- which(ring != 0)
    
    if (length(no_zero) >= 3) {
      ring_expr[i] <- sum(counts_matrix_row[index_bin]) / length(no_zero)
    }
  }
  
  return(ring_expr)
}


process_sub_matrix <- function(gene_group) {
  
  process_row <- function(counts_matrix_row) {
    # Filtra i geni con meno di tre counts nella matrice originale
    if (sum(counts_matrix_row != 0) >= 3) {
      # Calcola l'espressione media dei geni per ogni ring solo se soddisfa il filtro dei conteggi minimi
      gene_ring_expr_parallel(counts_matrix_row, rings)
    } else {
      rep(0, length(unique(rings)))  # Se non soddisfa il filtro, restituisci un vettore di zeri
    }
  }
  
  gene_ring_expr_values <- t(apply(counts_matrix[gene_group,], 1, process_row))
  return(gene_ring_expr_values)
}


# Funzione per calcolare l'espressione media dei geni per ogni ring in parallelo
calculate_gene_ring_expr_parallel <- function(gene_groups, counts_matrix, rings) {
  cluster <- makeCluster(num_cores)
  
  # Esporto la funzione e altre variabili necessarie nell'ambiente di esecuzione parallela
  
  
  clusterExport(cluster, c("gene_ring_expr_parallel", "counts_matrix", "rings", "process_sub_matrix"))
  gene_ring_expr_values <- parLapply(cluster, gene_groups, process_sub_matrix)
  stopCluster(cluster)
  
  return(gene_ring_expr_values)
}

# Calcolo l'espressione media dei geni per ogni ring in parallelo
gene_ring_expr_values <- calculate_gene_ring_expr_parallel(gene_groups, counts_matrix, rings)


# Costruisco la matrice matrice_rings
matrice_rings <- do.call(rbind, gene_ring_expr_values)
dimnames(matrice_rings) <- list(
  rownames(counts_matrix),
  paste0("ring", 1:ncol(matrice_rings))
  )  

# Rimuovo i geni che hanno meno di 5 rings con espressione diversa da zero
removed_genes <- apply(matrice_rings != 0, 1, sum) < 5

matrice_rings[matrice_rings==0] <- NA
correlations <- apply(matrice_rings, 1, cor, y=1:ncol(matrice_rings), method="spearman")
par(mfrow=c(3,3))
hist(correlations)
sorter <- order(abs(correlations), decreasing = T)
for(i in 1:8)
  plot(matrice_rings[sorter[i], ], xlab="Dist", ylab="Expr", main=rownames(matrice_rings)[i])

