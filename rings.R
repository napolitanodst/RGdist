
## requirements:
# BiocManager::install("SpatialExperiment")

config <- read.table("config.txt")[,2]
names(config) <- read.table("config.txt")[,1]

library(SpatialExperiment)
load(file.path(config["input_path"], "spe-spots.RData"))

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
num_cores <- 10
counts_matrix <- as.matrix(assays(spe)$counts)

# Divido gli indici dei geni in gruppi per parallelizzazione
gene_groups <- split(1:nrow(counts_matrix), rep(1:num_cores, each = ceiling(nrow(counts_matrix) / num_cores))[1:nrow(counts_matrix)])

# Funzione per calcolare l'espressione media dei geni per ogni ring
gene_ring_expr_parallel <- function(gene_id, counts_matrix, rings) {
  ring_expr <- numeric(length(unique(rings)))
  
  for (i in 1:length(unique(rings))) {
    index_bin <- which(unique(rings)[i] == rings)
    no_zero <- which(counts_matrix[gene_id, index_bin] != 0)
    
    if (length(no_zero) >= 3) {
      ring_expr[i] <- sum(counts_matrix[gene_id, index_bin]) / length(no_zero)
    } else {
      ring_expr[i] <- 0
    }
  }
  
  return(ring_expr)
}


# Funzione per calcolare l'espressione media dei geni per ogni ring in parallelo
calculate_gene_ring_expr_parallel <- function(gene_groups, counts_matrix, rings) {
  cluster <- makeCluster(num_cores)
  
  # Esporto la funzione e altre variabili necessarie nell'ambiente di esecuzione parallela
  clusterExport(cluster, c("gene_ring_expr_parallel", "counts_matrix", "rings"))
  
  gene_ring_expr_values <- parLapply(cluster, gene_groups, function(gene_indices) {
    gene_ring_expr_values <- lapply(gene_indices, function(gene_id) {
      # Filtra i geni con meno di tre counts nella matrice originale
      if (sum(counts_matrix[gene_id, ] != 0) >= 3) {
        # Calcola l'espressione media dei geni per ogni ring solo se soddisfa il filtro dei conteggi minimi
        gene_ring_expr_parallel(gene_id, counts_matrix, rings)
      } else {
        rep(0, length(unique(rings)))  # Se non soddisfa il filtro, restituisci un vettore di zeri
      }
    })
    return(gene_ring_expr_values)
  })
  stopCluster(cluster)
  
  # Unisci i risultati in un'unica lista
  gene_ring_expr_values <- unlist(gene_ring_expr_values)
  
  return(gene_ring_expr_values)
}

# Calcolo l'espressione media dei geni per ogni ring in parallelo
gene_ring_expr_values <- calculate_gene_ring_expr_parallel(gene_groups, counts_matrix, rings)


# Costruisco la matrice matrice_rings
matrice_rings <- matrix(unlist(gene_ring_expr_values), 
                        nrow = nrow(counts_matrix),
                        ncol = 15, 
                        byrow = TRUE,
                        dimnames = list(rownames(counts_matrix),
                                        c("ring1", "ring2", "ring3", "ring4", "ring5",
                                          "ring6", "ring7", "ring8", "ring9", "ring10",
                                          "ring11", "ring12", "ring13", "ring14", "ring15")))

# Rimuovo i geni che hanno meno di 5 rings con espressione diversa da zero
removed_genes <- apply(matrice_rings != 0, 1, sum) < 5


