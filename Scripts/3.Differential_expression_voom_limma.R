# ---------------------------------------
# 1) Paquetes
# ---------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma", "edgeR"), ask = FALSE)
install.packages("dplyr")

library(limma)
library(edgeR)
library(dplyr)

# ---------------------------------------
# 2) Rutas
# ---------------------------------------
input_path <- "C:/Users/victo/OneDrive/Documentos/AI4BIO_DATA/GSE279208/GSE279208_matriz_protein.txt"
output_path <- "C:/Users/victo/OneDrive/Documentos/AI4BIO_DATA/GSE279208/GSE279208_DE.txt"

# ---------------------------------------
# 3) Cargar matriz de expresión
# ---------------------------------------
expr <- read.table(
  input_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

genes <- expr$Symbol
expr_mat <- as.matrix(expr[, -1])
rownames(expr_mat) <- genes

# ---------------------------------------
# 4) Definir grupos (desde metadata)
# ---------------------------------------
control_samples <- c(
  "GSM8564259","GSM8564260","GSM8564261","GSM8564262","GSM8564263",
  "GSM8564264","GSM8564265","GSM8564266","GSM8564267","GSM8564268"
)

infection_samples <- c(
  "GSM8564269","GSM8564270","GSM8564271","GSM8564272","GSM8564273",
  "GSM8564274","GSM8564275","GSM8564276","GSM8564277","GSM8564278",
  "GSM8564279","GSM8564280","GSM8564281","GSM8564282","GSM8564283",
  "GSM8564284"
)

group <- factor(
  ifelse(colnames(expr_mat) %in% control_samples, "Control", "Infection"),
  levels = c("Control", "Infection")
)

# Comprobación rápida
table(group)

# ---------------------------------------
# 5) edgeR + voom
# ---------------------------------------
dge <- DGEList(counts = expr_mat, group = group)

# Filtrado mínimo recomendado
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Infection_vs_Control")

v <- voom(dge, design, plot = TRUE)

# ---------------------------------------
# 6) limma
# ---------------------------------------
fit <- lmFit(v, design)
fit <- eBayes(fit)

deg <- topTable(
  fit,
  coef = "Infection_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

deg$Symbol <- rownames(deg)

# ---------------------------------------
# 7) Guardar resultados
# ---------------------------------------
write.table(
  deg,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Genes analizados:", nrow(deg), "\n")
cat("Archivo guardado en:\n", output_path, "\n")
