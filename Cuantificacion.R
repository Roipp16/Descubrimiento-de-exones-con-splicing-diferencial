files <- list.files("/home/roi/Downloads/gene_counts/", pattern = "\\.sf$", full.names = TRUE)
gene_counts_list <- lapply(files, function(x) read.delim(x, sep = "\t"))
sample_names <- gsub("_quant.genes", "", tools::file_path_sans_ext(basename(files)))

library(dplyr)
library(purrr)

gene_counts_list <- lapply(seq_along(files), function(i) {
  df <- read.delim(files[i], sep = "\t")
  df <- df[, c("Name", "NumReads")]  
  colnames(df)[2] <- sample_names[i]  
  return(df)
})

counts_by_sample <- purrr::reduce(gene_counts_list, dplyr::full_join, by = "Name")
rownames(counts_by_sample) <- counts_by_sample$Name
rm(gene_counts_list)
rm(files)
rm(sample_names)
gc()

library(edgeR)
DGEList<- DGEList(counts_by_sample, remove.zeros = TRUE)

groups <- c("interm", "high", "interm", "low", "interm", "high", "low", "low", "high")
clinical_status <- data.frame(
  sample = colnames(counts_by_sample)[-1], 
  group = factor(groups, levels = c("low", "interm", "high"))
)
design <- model.matrix(~0 + group, data = clinical_status)
colnames(design) <- gsub("group", "", colnames(design))

keep.expr<-filterByExpr(DGEList,design)
dge<-DGEList[keep.expr,]

log_counts <- log2(dge$counts + 1)

dge<- calcNormFactors(dge, method = "TMM")
tmmexp<-cpm(dge, log = T, prior.count = 3)

hist(tmmexp, 100)
boxplot(tmmexp)
boxplot(log_counts)


rm(DGEList)
rm(keep.expr)
rm(counts_by_sample)


fit<-lmFit(tmmexp,design)

cont_mat <- makeContrasts(
  high_vs_low = high - low,

  levels = design
)

fit2_HL <- contrasts.fit(fit, cont_mat)
fit2_HL <- eBayes(fit2_HL)
plotSA(fit2_HL)


de_res_HL <- topTable(fit2_HL, number = Inf)
de_res_HL$gene_symbol <- rownames(de_res_HL)
signif_genes_HL <- de_res_HL$gene_symbol[de_res_HL$P.Value < 0.05]
write.csv(de_res_HL, file = "high_vs_low", row.names = TRUE)

dt_HL<-decideTests(fit2_HL, adjust.method = "none", p.value = 0.05)


fit<-lmFit(tmmexp,design)
cont_mat <- makeContrasts(
  high_vs_interm = high - interm,
  levels = design
)

fit2_HI <- contrasts.fit(fit, cont_mat)
fit2_HI <- eBayes(fit2_HI)
plotSA(fit2_HI)
de_res_HI <- topTable(fit2_HI, number = Inf)
de_res_HI$gene_symbol <- rownames(de_res_HI)
signif_genes_HI <- de_res_HI$gene_symbol[de_res_HI$P.Value < 0.05]
write.csv(de_res_HI, file = "high_vs_interm", row.names = TRUE)




fit<-lmFit(tmmexp,design)
cont_mat <- makeContrasts(
  interm_vs_low = interm - low,
  levels = design
)
fit2_IL <- contrasts.fit(fit, cont_mat)
fit2_IL <- eBayes(fit2_IL)
plotSA(fit2_IL)
de_res_IL <- topTable(fit2_IL, number = Inf)

de_res_IL$gene_symbol <- rownames(de_res_IL)
signif_genes_IL <- de_res_IL$gene_symbol[de_res_IL$P.Value < 0.05]
write.csv(de_res_IL, file = "interm_vs_low", row.names = TRUE)

common_genes_sig2 <-intersect(intersect(signif_genes_HI, signif_genes_HL), signif_genes_IL)


dt_IL<- decideTests(fit2_IL, adjust.method = "none", p.value = 0.05)
dt_HI<- decideTests(fit2_HI, adjust.method = "none", p.value = 0.05)
dt_HL<-decideTests(fit2_HL, adjust.method = "none", p.value = 0.05)
  
library(ComplexHeatmap)
sig <- list(
  HL = signif_genes_HL,
  HI = signif_genes_HI,
  IL = signif_genes_IL
)
mat <- list_to_matrix(sig)
library(ComplexUpset)
library(ggplot2)
df <- as.data.frame(mat)
df$gene <- rownames(mat)


upset(
  df,
  intersect = c("HL", "HI", "IL"),
  name = "Genes",
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 3)
    ) + 
      labs(y = "Número de genes")
  ),
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(
        panel.grid = element_blank()
      ),
      'Intersection size' = theme(
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1)
      )
    )
  ),
  width_ratio = 0.2,
  min_size = 1,
  sort_intersections_by = 'degree'
)


vennDiagram(
  cbind(dt_HI, dt_HL, dt_IL),include=c("both") ,
  counts.col = c("red", "blue"),
  circle.col = c("#FF0000", "#00B200", "#0000FF"),   
  circle.lwd = 2,
  cex = 1.2,
  cat.cex = 1.5,
  cat.col = c("#FF0000", "#00B200", "#0000FF"),
  names = c("High vs Low", "High vs Interm", "Interm vs Low"),
  cat.pos = c(-15, 15, 135),
  cat.dist = 0.07,
  main = "Venn Diagram",
  main.cex = 1.6,
)

vennDiagram(
  cbind(dt_HI, dt_HL, dt_IL),include=c("up", "down") ,
  counts.col = c("red", "blue"),
  circle.col = c("#FF0000", "#00B200", "#0000FF"),   
  circle.lwd = 2,
  cex = 1.2,
  cat.cex = 1.5,
  cat.col = c("#FF0000", "#00B200", "#0000FF"),
  names = c("High vs Low", "High vs Interm", "Interm vs Low"),
  cat.pos = c(-15, 15, 135),
  cat.dist = 0.07,
  main = "Venn Diagram",
  main.cex = 1.6,
)

#volcano plot 1
de_res_HL$diffexpressed <- "NO"
de_res_HL$diffexpressed[de_res_HL$logFC > 0.6 & de_res_HL$P.Value < 0.05] <- "UP"
de_res_HL$diffexpressed[de_res_HL$logFC < -0.6 & de_res_HL$P.Value < 0.05] <- "DOWN"

de_res_HL$delabel <- NA
de_res_HL$delabel[de_res_HL$diffexpressed != "NO"] <- de_res_HL$gene_symbol[de_res_HL$diffexpressed != "NO"]

mycolors <- c("DOWN" = "blue", "UP" = "red", "NO" = "black")

library(ggplot2)
library(ggrepel)

ggplot(de_res_HL, aes(x = logFC, y = -log10(P.Value), color = diffexpressed, label = delabel)) +
  geom_point() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  scale_color_manual(values = mycolors) +
  geom_text_repel(size = 3, max.overlaps = 10) +
  theme_minimal() + labs(title = "High vs Low")



#volcano plot 2
de_res_HI$diffexpressed <- "NO"
de_res_HI$diffexpressed[de_res_HI$logFC > 0.6 & de_res_HI$P.Value < 0.05] <- "UP"
de_res_HI$diffexpressed[de_res_HI$logFC < -0.6 & de_res_HI$P.Value < 0.05] <- "DOWN"

de_res_HI$delabel <- NA
de_res_HI$delabel[de_res_HI$diffexpressed != "NO"] <- de_res_HI$gene_symbol[de_res_HI$diffexpressed != "NO"]

mycolors <- c("DOWN" = "blue", "UP" = "red", "NO" = "black")

library(ggplot2)
library(ggrepel)

ggplot(de_res_HI, aes(x = logFC, y = -log10(P.Value), color = diffexpressed, label = delabel)) +
  geom_point() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  scale_color_manual(values = mycolors) +
  geom_text_repel(size = 3, max.overlaps = 10) +
  theme_minimal() + labs(title = "High vs Interm")


#volcano plot 3

de_res_IL$diffexpressed <- "NO"
de_res_IL$diffexpressed[de_res_IL$logFC > 0.6 & de_res_IL$P.Value < 0.05] <- "UP"
de_res_IL$diffexpressed[de_res_IL$logFC < -0.6 & de_res_IL$P.Value < 0.05] <- "DOWN"

de_res_IL$delabel <- NA
de_res_IL$delabel[de_res_IL$diffexpressed != "NO"] <- de_res_IL$gene_symbol[de_res_IL$diffexpressed != "NO"]

mycolors <- c("DOWN" = "blue", "UP" = "red", "NO" = "black")

library(ggplot2)
library(ggrepel)

ggplot(de_res_IL, aes(x = logFC, y = -log10(P.Value), color = diffexpressed, label = delabel)) +
  geom_point() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  scale_color_manual(values = mycolors) +
  geom_text_repel(size = 3, max.overlaps = 10) +
  theme_minimal()

#creamos df para con significativos y p valores para los tres grupos y luego juntamos
sigHI <- de_res_HI[de_res_HI$P.Value < 0.05, ]
sigHI <- data.frame(Gene = rownames(sigHI), P.Value = sigHI$P.Value)

sigHL <- de_res_HL[de_res_HL$P.Value < 0.05, ]
sigHL <- data.frame(Gene = rownames(sigHL), P.Value = sigHL$P.Value)

sigIL <- de_res_IL[de_res_IL$P.Value < 0.05, ]
sigIL <- data.frame(Gene = rownames(sigIL), P.Value = sigIL$P.Value)

hi <- setNames(-log10(sigHI$P.Value), sigHI$Gene)
hl <- setNames(-log10(sigHL$P.Value), sigHL$Gene)
il <- setNames(-log10(sigIL$P.Value), sigIL$Gene)

all_genes <- unique(c(names(hi), names(hl), names(il)))

heatmap_mat <- data.frame(
  HI = hi[all_genes],
  HL = hl[all_genes],
  IL = il[all_genes]
)

heatmap_mat[is.na(heatmap_mat)] <- 0
rownames(heatmap_mat) <- all_genes

library(ComplexHeatmap)
library(circlize)

min_p <- apply(heatmap_mat, 1, max)  # mayor -log10(p) = menor p
top_genes <- names(sort(min_p, decreasing = TRUE))
heatmap_subset <- heatmap_mat[top_genes, ]

# Convertir a matriz
heatmap_matrix <- as.matrix(heatmap_subset)

# Colores
col_fun <- colorRamp2(c(0, max(heatmap_matrix)), c("white", "red"))

Heatmap(heatmap_matrix,
        name = "-log10(p)",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_dend_width = unit(2, "cm"),    # ancho dendrograma filas (por defecto es menor)
        column_dend_height = unit(1, "cm"))



genes <- c("UGCG", "B4GALT5", "B4GALT6", "ST3GAL5", "ST8SIA1", "ST8SIA3", "ST8SIA5", "B4GALT1", "B3GALT1", "B3GALT4", "ST3GAL2", "ST3GAL3")

sigHI <- de_res_HI[de_res_HI$P.Value < 0.05 & rownames(de_res_HI) %in% genes, ]
sigHI <- data.frame(Gene = rownames(sigHI), P.Value = sigHI$P.Value)

sigHL <- de_res_HL[de_res_HL$P.Value < 0.05 & rownames(de_res_HL) %in% genes, ]
sigHL <- data.frame(Gene = rownames(sigHL), P.Value = sigHL$P.Value)

sigIL <- de_res_IL[de_res_IL$P.Value < 0.05 & rownames(de_res_IL) %in% genes, ]
sigIL <- data.frame(Gene = rownames(sigIL), P.Value = sigIL$P.Value)

selec1 <- setNames(sigHI$P.Value, sigHI$Gene)
selec2 <- setNames(sigHL$P.Value, sigHL$Gene)
selec3 <- setNames(sigIL$P.Value, sigIL$Gene)

all_selec <- unique(c(names(selec1), names(selec2), names(selec3)))

df_selec <- data.frame(
  HI = selec1[all_selec],
  HL = selec2[all_selec],
  IL = selec3[all_selec]
)

df_selec[is.na(df_selec)] <- 1  # relleno con 1 (p no significativo)

df_selec <- df_selec[genes, , drop = FALSE] #solo hay 2 genes?


conteo<- genes%in%rownames(counts_by_sample)#al inicio estan todos menos B4GALTN1, pero B4GALT1 sí está
"B4GALTN1"%in% rownames(counts_by_sample)#FALSE, no está
"B4GALT1"%in% rownames(counts_by_sample)#TRUE


#vamos a buscar donde se pierden: 
genes%in%signif_genes_HI | genes%in%signif_genes_HL | genes%in%signif_genes_IL #solo hay tres de significativos
c("B4GALT6", "ST8SIA3", "B4GALT1") %in%signif_genes_HI | genes%in%signif_genes_HL | genes%in%signif_genes_IL #son estos 
gogenes <-c("B4GALT6", "ST8SIA3", "B4GALT1")


#vamos a enriquecer
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(clusterProfiler)
library(ggplot2)
library(DOSE)
entrez_genes <- bitr(gogenes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene         = entrez_genes$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "ALL",         
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

    
df_ego <- ego@result
ggplot(df_ego, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "GO for B4GALT6, ST8SIA3, B4GALT1", x = "GO Term", y = "Gene Count") +
  theme_minimal()

KEGG_genes <- enrichKEGG(gene         = entrez_genes$ENTREZID, organism     = "hsa", pvalueCutoff = 0.05)
top_kegg<- KEGG_genes@result

ggplot(top_kegg, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "KEGG for B4GALT6, ST8SIA3, B4GALT1",
       x = "KEGG Pathway",
       y = "Gene Count",
       fill = "p-value") +
  theme_minimal()



