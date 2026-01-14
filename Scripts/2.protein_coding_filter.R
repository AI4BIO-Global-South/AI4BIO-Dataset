# Instalar Bioconductor (si no lo tienes)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar biomaRt
BiocManager::install("biomaRt")

# Instalar dplyr
install.packages("dplyr")
library(biomaRt)
library(dplyr)

# ---------------------------------------
# 1) Definir rutas
# ---------------------------------------
input_path <- "C:/Users/victo/OneDrive/Documentos/AI4BIO_DATA/GSE25226/GSE25226_matriz.txt"

output_path <- "C:/Users/victo/OneDrive/Documentos/AI4BIO_DATA/GSE25226/GSE25226_matriz_protein.txt"

# ---------------------------------------
# 2) Cargar matriz de expresión
#    Primera columna: Symbol
# ---------------------------------------
expr <- read.table(
  input_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ---------------------------------------
# 3) Eliminar genes sin símbolo
# ---------------------------------------
expr <- expr %>%
  filter(!is.na(Symbol) & Symbol != "")

# ---------------------------------------
# 4) Conectar a Ensembl (biomaRt)
# ---------------------------------------
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

annotation <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters = "hgnc_symbol",
  values = unique(expr$Symbol),
  mart = mart
)

# ---------------------------------------
# 5) Obtener genes protein-coding
# ---------------------------------------
protein_coding_genes <- annotation %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(hgnc_symbol) %>%
  unique()

# ---------------------------------------
# 6) Filtrar matriz de expresión
# ---------------------------------------
expr_protein <- expr %>%
  filter(Symbol %in% protein_coding_genes)

# ---------------------------------------
# 7) Guardar matriz filtrada
# ---------------------------------------
write.table(
  expr_protein,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---------------------------------------
# 8) Resumen
# ---------------------------------------
cat("Genes totales:", nrow(expr), "\n")
cat("Genes protein-coding:", nrow(expr_protein), "\n")
cat("Archivo guardado en:\n", output_path, "\n")

