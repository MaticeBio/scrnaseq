#!/usr/bin/env Rscript

# to use nf variables: "${meta.id}"

# load libraries
library(anndataR)
library(SeuratObject)
library(SingleCellExperiment)

# read input
adata <- read_h5ad("${h5ad}")

# Seurat (SeuratObject) mangles feature names containing '_' to '-' at object
# creation, after which anndataR re-keys the counts layer by the ORIGINAL
# (underscore) var_names -> "No feature overlap between existing object and new
# layer data" and the conversion aborts. Pre-sanitize var_names so anndataR and
# Seurat agree (and Seurat/SCE rownames stay identical). This bites references
# whose gene IDs contain '_' (e.g. C. elegans RefSeq CELE_...); it is a no-op for
# Ensembl IDs (ENSG.../ENSMUSG..., which never contain '_').
adata\$var_names <- gsub("_", "-", adata\$var_names)

# convert to Seurat
obj <- adata\$as_Seurat(x_mapping = "counts")

# save files
dir.create(file.path("$meta.id"), showWarnings = FALSE)
saveRDS(obj, file = "${meta.id}_${meta.input_type}_matrix.seurat.rds")

# convert to SingleCellExperiment. x_mapping="counts" so the primary assay is
# named "counts" (the default leaves it as X/count and breaks assay(obj,"counts")).
obj <- adata\$as_SingleCellExperiment(x_mapping = "counts")

# save files
dir.create(file.path("$meta.id"), showWarnings = FALSE)
saveRDS(obj, file = "${meta.id}_${meta.input_type}_matrix.sce.rds")

#
# save versions file
#
versions_file <- file("versions.yml")
write(
    paste(
        '${task.process}:',
        paste0('  r-base: "', R.Version()\$version.string, '"'),
        paste0('  anndataR: "', as.character(packageVersion("anndataR")), '"'),
        paste0('  SeuratObject: "', as.character(packageVersion("SeuratObject")), '"'),
        paste0('  SingleCellExperiment: "', as.character(packageVersion("SingleCellExperiment")), '"'),
        sep = "\\n"
    ),
    versions_file
)
close(versions_file)
