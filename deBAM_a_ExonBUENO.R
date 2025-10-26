### Carga de paquetes necesarios de R en el servidor
lib.path="/gpfs/res_projects/fsjd3/rlib/"
.libPaths(c(lib.path,.libPaths()))
R_LIBS_USER="/gpfs/res_projects/fsjd3/rlib/"

### Para este script, se usan diferentes checkpoints ya que es un script computacionalmente exigente
### Como inputs es necesario disponer de los ficheros BAM de los datos a analizar.

library(GenomicFeatures) 
library(txdbmaker)
library(DEXSeq)
library(DESeq2)
library(Rsamtools)
library(GenomicAlignments)
library(GenomeInfoDb)
library(BiocParallel)

# Nombres de archivos por si hay algun paso intermedio ya procesado
rdata_file <- "bam.Rdata"
rdata_file2 <- "diff_exons_finished.Rdata"


# 1. Carga BAM files, anotación y cálculo de Overlaps. El input serán los BAM y un archivo genes.gtf
### Se consulta si tenemos ya los BAM analizados: 

if (file.exists(rdata_file)) {
  cat("Cargando datos desde", rdata_file, "...\n")
  load(rdata_file, envir = parent.frame(), verbose = FALSE)
} else {
  cat("No existe", rdata_file, "- generando...\n")
  
  txdb <- makeTxDbFromGFF("../genes.gtf")
  flattenedAnnotation <- exonicParts(txdb, linked.to.single.gene.only=TRUE)
  names(flattenedAnnotation) <- sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)

  ruta_bam <- "../NBtest/FINAL/star_salmon/"
  mis_bam <- list.files(path = ruta_bam, pattern = "\\.bam$", full.names = TRUE)
  bamFiles <- BamFileList(mis_bam)

  bam_seqinfo <- seqinfo(bamFiles[[1]])
  seqlevelsStyle(flattenedAnnotation) <- seqlevelsStyle(bam_seqinfo)
  bam_seqlevels <- seqlevels(bam_seqinfo)
  seqlevels(flattenedAnnotation, pruning.mode = "coarse") <- bam_seqlevels
  seqinfo(flattenedAnnotation) <- bam_seqinfo

  register(SerialParam())
  
  ### Para optimizar el tiempo de ejecución se separan los BAMs de 3 en 3, ya que en mayor número
  ### se requiriría demasiado poder computacional. Con 3 no hay problema
  
  batch_size <- 3
  n_bams <- length(bamFiles)
  batch_indices <- split(seq_len(n_bams), ceiling(seq_len(n_bams) / batch_size))
  se_list <- list()

  for (i in seq_along(batch_indices)) {
    cat("Procesando lote", i, "...\n")
    idx <- batch_indices[[i]]
    bam_batch <- bamFiles[idx]

    se_batch <- summarizeOverlaps(
      flattenedAnnotation,
      bam_batch,
      singleEnd = FALSE,
      fragments = TRUE,
      ignore.strand = TRUE
    )

    se_list[[i]] <- se_batch
  }
  se_combined <- do.call(cbind, se_list)
  
  save(se_combined, file = rdata_file)
}
# -------------      FIN DEL CHECKPOINT 1      -------------


# 2. Cálculo DEXSeq dataset y análisis comparativo
### Se comprueba existencia del checkpoint 2:

if (file.exists(rdata_file2)) {
  cat("Cargando datos desde", rdata_file2, "...\n")
  load(rdata_file2, envir = parent.frame(), verbose = FALSE)
} else {
  cat("No existe", rdata_file2, "- generando análisis DEXSeq...\n")
  
  groups <- c("interm", "high", "interm", "low", "interm", "high", "low", "low", "high")
  colData(se_combined)$sample <- factor(colnames(se_combined))
  
### IMPORTANTE: el que haya en primera posición actuará como referencia. hacemos 1 para cada grupo
  #Filtramos para hacer las comparaciones: 
  
  # 1r grupo low-high:
  cat("Empezando comparación low x high\n")
  idx <- which(groups %in% c("low", "high"))
  se_sub <- se_combined[, idx]
  groups_sub <- groups[idx]
  
  colData(se_sub)$sample <- factor(colnames(se_sub))
  colData(se_sub)$condition <- factor(groups_sub, levels = c("low", "high"))
  
  dxd_low_high <- DEXSeqDataSetFromSE(se_sub, design = ~ sample + exon + condition:exon)
  dxd_low_high <- estimateSizeFactors(dxd_low_high)
  dxd_low_high <- estimateDispersions(dxd_low_high)
  dxd_low_high <- testForDEU(dxd_low_high)
  dxd_low_high <- estimateExonFoldChanges(dxd_low_high, fitExpToVar = "condition")

  ### como hablamos de mucho tiempo de ejecución, se guarda un subcheckpoint: 
  save.image(file="1diff_exons.RData")
  cat("Rdata saved in 1diffexons")
  
  # resultados
  res_low_high <- DEXSeq::DEXSeqResults(dxd_low_high)
  res_low_high_df <- as.data.frame(res_low_high)
  
  # Fuerza a que toda columna lista se vuelva carácter para evitar errores
  for (col in names(res_low_high_df)) {
    if (is.list(res_low_high_df[[col]])) {
      res_low_high_df[[col]] <- vapply(res_low_high_df[[col]], 
                                       function(x) paste(as.character(x), collapse = ";"), 
                                       character(1))
    }
  }
  
  # Exportar a CSV
  write.csv(res_low_high_df, file = "DEXSEQ_low_vs_high.csv", row.names = TRUE)
  cat("csv comparacion LOW-HIGH guardado\n")



  
  # 2º grupo interm-low
  cat("Empezando comparación interm x low\n")
  idx_interm_low <- which(groups %in% c("interm", "low"))
  se_interm_low <- se_combined[, idx_interm_low]
  groups_interm_low <- groups[idx_interm_low]
  
  colData(se_interm_low)$sample <- factor(colnames(se_interm_low))
  colData(se_interm_low)$condition <- factor(groups_interm_low, levels = c("interm", "low"))
  
  dxd_interm_low <- DEXSeqDataSetFromSE(se_interm_low, design = ~ sample + exon + condition:exon)
  dxd_interm_low <- estimateSizeFactors(dxd_interm_low)
  dxd_interm_low <- estimateDispersions(dxd_interm_low)
  dxd_interm_low <- testForDEU(dxd_interm_low)
  dxd_interm_low <- estimateExonFoldChanges(dxd_interm_low, fitExpToVar = "condition")
   
  # resultados
  res_interm_low <- DEXSeq::DEXSeqResults(dxd_interm_low)
  res_interm_low_df <- as.data.frame(res_interm_low)
  
  # Fuerza a que toda columna lista se vuelva carácter
  for (col in names(res_interm_low_df)) {
    if (is.list(res_interm_low_df[[col]])) {
      res_interm_low_df[[col]] <- vapply(res_interm_low_df[[col]], 
                                         function(x) paste(as.character(x), collapse = ";"), 
                                         character(1))
    }
  }
  
  # Exportar a CSV
  write.csv(res_interm_low_df, file = "DEXSEQ_interm_vs_low.csv", row.names = TRUE)
  cat("csv comparacion INTERM-LOW guardado\n")

  
  # 3er grupo high-interm
  cat("Empezando comparación high x interm\n")
  idx_high_interm <- which(groups %in% c("high", "interm"))
  se_high_interm <- se_combined[, idx_high_interm]
  groups_high_interm <- groups[idx_high_interm]
  
  colData(se_high_interm)$sample <- factor(colnames(se_high_interm))
  colData(se_high_interm)$condition <- factor(groups_high_interm, levels = c("high", "interm"))
  
  dxd_high_interm <- DEXSeqDataSetFromSE(se_high_interm, design = ~ sample + exon + condition:exon)
  dxd_high_interm <- estimateSizeFactors(dxd_high_interm)
  dxd_high_interm <- estimateDispersions(dxd_high_interm)
  dxd_high_interm <- testForDEU(dxd_high_interm)
  dxd_high_interm <- estimateExonFoldChanges(dxd_high_interm, fitExpToVar = "condition")
   
  # resultados
  res_high_interm <- DEXSeq::DEXSeqResults(dxd_high_interm)
  res_high_interm_df <- as.data.frame(res_high_interm)
  
  # Fuerza a que toda columna lista se vuelva carácter
  for (col in names(res_high_interm_df)) {
    if (is.list(res_high_interm_df[[col]])) {
      res_high_interm_df[[col]] <- vapply(res_high_interm_df[[col]], 
                                          function(x) paste(as.character(x), collapse = ";"), 
                                          character(1))
    }
  }
  
  # Exportar a CSV
  write.csv(res_high_interm_df, file = "DEXSEQ_high_vs_interm.csv", row.names = TRUE)
  cat("csv comparacion HIGH-INTERM guardado\n")
  
  # guardar objetos y sesión
  save(dxd_low_high, dxd_interm_low, dxd_high_interm, file = rdata_file2)
  save.image(file = "diff_exons_finished.RData")
  
}
