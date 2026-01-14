library(dplyr)

deg <- read.table(
  "C:/Users/victo/OneDrive/Documentos/AI4BIO_DATA/GSE51808/GSE51808_DE.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ===============================
# 1. Eliminar columna Entrez (mal llamada Symbol)
# ===============================

deg <- deg %>%
  select(-Symbol)

# ===============================
# 2. Renombrar ID -> Symbol
# ===============================

deg <- deg %>%
  rename(Symbol = ID)

# ===============================
# 3. Reordenar columnas (opcional pero recomendado)
# ===============================

deg <- deg %>%
  select(logFC, AveExpr, t, P.Value, adj.P.Val, B, Symbol)

# ===============================
# 4. Colapsar genes repetidos
# ===============================

deg_collapsed <- deg %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%
  group_by(Symbol) %>%
  slice(1) %>%
  ungroup()

# ===============================
# 5. Guardar archivo final
# ===============================

write.table(
  deg_collapsed,
  file = "C:/Users/victo/OneDrive/Documentos/AI4BIO_DATA/GSE51808/GSE51808_DE_clean.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ===============================
# 6. Chequeos
# ===============================

any(duplicated(deg_collapsed$Symbol))  # debe ser FALSE
head(deg_collapsed$Symbol)



