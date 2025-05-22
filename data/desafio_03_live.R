
library(tidyverse)
library(pheatmap)

ensembl_to_genename <- read_tsv(file = "./data/ensembl_genes.tsv")
oncogene_names <- c(
  "EGFR",
  "TP53",
  "MYC",
  "KRAS",
  "BRAF",
  "PIK3CA",
  "CDKN2A",
  "MDM2",
  "RB1",
  "ERBB2"
)

oncogene_ensembl_id <- ensembl_to_genename %>%
  filter(
    hgnc_symbol %in% oncogene_names
  )


# Callback para filtrar os chunks (EGFR == ENSG00000146648)
callback <- function(x, pos) {
  x %>% filter(Gene %in% oncogene_ensembl_id$ensembl_gene_id)
}

# Leitura em chunks
result <- read_tsv_chunked(
  file = "/home/oandrefonseca/Downloads/rna_cancer_sample.tsv",
  callback = DataFrameCallback$new(callback),
  chunk_size = 5000 # Reduzir em caso de pouca memoria
)

# Combina resultados dos chunks (caso múltiplos batam a condição)
final_result <- bind_rows(result)

write_tsv(final_result, file = "./data/oncogene_heatmap_data.tsv")

#

result <- read_tsv(file = "./data/oncogene_heatmap_data.tsv")
result <- result %>%
  mutate(
    BOLINHA = stringr::str_remove_all(Cancer, "\\s+\\(\\w+\\)")
  )

#\\s+ = um ou mais espacos em branco
#\\( = escape e encontre um (
#\\w+ = w+ um ou mais caracteres
#\\) = escape e encontre um )

#

sample_annotation <- result %>%
  select(Sample, BOLINHA) %>%
  distinct() %>%
  column_to_rownames(var = "Sample")

gene_avg_tumor <- result %>%
  group_by(Gene, BOLINHA) %>%
  summarise(
    avg_pTPM = mean(pTPM, trim = 0.1)
  )

#

gene_sample <- result %>%
  select(Gene, Sample, pTPM) %>%
  pivot_wider(
    names_from = "Sample",
    values_from = "pTPM"
  )

gene_sample <- gene_sample %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

#

pheatmap(
  gene_sample,
  scale = "column",
  cluster_cols = FALSE,
  annotation_col = sample_annotation,
  show_colnames = FALSE
  )

#

gene_sample_transpose <- gene_sample %>%
  t()

# ABORTADO

pheatmap(
  gene_sample_transpose,
  #annotation_row = sample_annotation,
  kmeans_k = 10
)

#

gene_avg_tumor <- gene_avg_tumor %>%
  select(Gene, Cancer, avg_pTPM) %>%
  pivot_wider(
    names_from = "Cancer",
    values_from = "avg_pTPM"
  )

gene_avg_tumor <- gene_avg_tumor %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

#

pheatmap(
  gene_avg_tumor,
  scale = "row"
)
