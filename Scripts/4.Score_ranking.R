# ---------------------------------------
# 1) Paquetes
# ---------------------------------------
install.packages("dplyr")  # solo si no está instalado
library(dplyr)

# ---------------------------------------
# 2) Rutas
# ---------------------------------------
input_path <- "C:/Users/victo/OneDrive/Documentos/AI4BIO_DATA/GSE40628/GSE40628_DE_clean.txt"

output_ranked <- "C:/Users/victo/OneDrive/Documentos/AI4BIO_DATA/GSE40628/GSE40628_DE2_ranking_score.txt"
output_top1000 <- "C:/Users/victo/OneDrive/Documentos/AI4BIO_DATA/GSE40628/GSE40628_DE2_top1000.txt"

# ---------------------------------------
# 3) Cargar tabla DE
# ---------------------------------------
deg <- read.table(
  input_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# ---------------------------------------
# 4) Verificación mínima
# ---------------------------------------
required_cols <- c("logFC", "P.Value", "Symbol")
if (!all(required_cols %in% colnames(deg))) {
  stop("El archivo no contiene las columnas necesarias: logFC, P.Value, Symbol")
}

# ---------------------------------------
# 5) Calcular score
#    score = |logFC| * -log10(P.Value)
# ---------------------------------------
deg <- deg %>%
  mutate(score = abs(logFC) * (-log10(P.Value)))

# ---------------------------------------
# 6) Ordenar TODOS los genes por score
# ---------------------------------------
deg_ranked <- deg %>%
  arrange(desc(score))

# ---------------------------------------
# 7) Top 1000 (subset del ranking)
# ---------------------------------------
top1000 <- deg_ranked %>%
  slice_head(n = 1000)

# ---------------------------------------
# 8) Guardar resultados
# ---------------------------------------
write.table(
  deg_ranked,
  file = output_ranked,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  top1000,
  file = output_top1000,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---------------------------------------
# 9) Resumen
# ---------------------------------------
cat("Genes totales:", nrow(deg_ranked), "\n")
cat("Top 1000:", nrow(top1000), "\n")
cat("Ranking completo guardado en:\n", output_ranked, "\n")
cat("Top 1000 guardado en:\n", output_top1000, "\n")

