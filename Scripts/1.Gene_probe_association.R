############################################################
## ANOTAR MATRIZ DE EXPRESIÓN CON GENE SYMBOL (GPL2700)
## GSE28405 — SCRIPT DEFINITIVO
############################################################

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

## Rutas
expr_file  <- "C:/Users/victo/OneDrive/Desktop/Results-AI4BIO/GSE28405/GSE28405_subset.txt"
annot_file <- "C:/Users/victo/OneDrive/Desktop/Results-AI4BIO/GSE28405/GPL2700.annot/GPL2700.annot"
output_file <- "C:/Users/victo/OneDrive/Desktop/Results-AI4BIO/GSE28405/GSE28405_expression_annotated.txt"

############################################################
## 1. LEER MATRIZ DE EXPRESIÓN
############################################################

expr <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## Verificación
if (!"ID_REF" %in% colnames(expr)) {
  stop("ERROR: La matriz de expresión no contiene la columna ID_REF")
}

############################################################
## 2. LEER ANOTACIÓN GPL2700 (FORMATO GEO)
############################################################

lines <- readLines(annot_file)

start <- grep("^!platform_table_begin", lines)
end   <- grep("^!platform_table_end", lines)

if (length(start) == 0 | length(end) == 0) {
  stop("ERROR: No se encontró la tabla en GPL2700.annot")
}

table_lines <- lines[(start + 1):(end - 1)]

annot <- read.table(
  text = table_lines,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  check.names = FALSE
)

if (!all(c("ID", "Gene symbol") %in% colnames(annot))) {
  stop("ERROR: La anotación no contiene ID y Gene symbol")
}

############################################################
## 3. SUBSET DE ANOTACIÓN (SOLO LO NECESARIO)
############################################################

annot_sub <- annot %>%
  select(
    ID_REF = ID,
    Symbol = `Gene symbol`
  )

############################################################
## 4. JOIN
############################################################

expr_annotated <- expr %>%
  left_join(annot_sub, by = "ID_REF") %>%
  relocate(Symbol, ID_REF)

############################################################
## 5. RESUMEN
############################################################

cat("Total de sondas:", nrow(expr_annotated), "\n")
cat("Con Symbol:", sum(!is.na(expr_annotated$Symbol)), "\n")
cat("Sin Symbol:", sum(is.na(expr_annotated$Symbol)), "\n")

############################################################
## 6. GUARDAR RESULTADO
############################################################

write.table(
  expr_annotated,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Archivo final guardado en:\n", output_file, "\n")
