library(ComplexUpset)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(limma)
library(ggrepel)
library(tibble)
library(edgeR)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

#setwd("/home/roi/Downloads/bam/analisis/analisis/")
#Partimos de 3 tablas normalizadas de conteo de exones
LH_df<- read.csv("./DEXSEQ_low_vs_high.csv", sep = ",")
IL_df<- read.csv("./DEXSEQ_interm_vs_low.csv", sep = ",")
HI_df<- read.csv("./DEXSEQ_high_vs_interm.csv", sep = ",")

#TABLA COMP LH: extraemos conteos x muestra
count_cols <- grep("^countData", colnames(LH_df))
count_mat <- LH_df[, count_cols]
rownames(count_mat) <- LH_df$X

#empezamos el analisis con la nueva matriz
DGEList<- DGEList(count_mat, remove.zeros = TRUE)
groups <- c("high", "low","high", "low", "low", "high")
clinical_status <- data.frame(
  sample = colnames(count_mat), 
  group = factor(groups, levels = c("low", "interm", "high"))
)
design <- model.matrix(~0 + group, data = clinical_status)

keep.expr<-filterByExpr(DGEList,design)
dge<-DGEList[keep.expr,]

log_counts <- log2(dge$counts + 1)

dge<- calcNormFactors(dge, method = "TMM")
tmmexp<-cpm(dge, log = T, prior.count = 3)

fit<-lmFit(tmmexp,design)

cont_mat <- makeContrasts(
  high_vs_low = grouphigh - grouplow,
  
  levels = design
)

fit2_HL <- contrasts.fit(fit, cont_mat)
fit2_HL <- eBayes(fit2_HL)

de_res_HL <- topTable(fit2_HL, number = Inf)
de_res_HL$ENSEMBL <- rownames(de_res_HL)

de_res_HL$ENSEMBL<- stringr::str_split_fixed(de_res_HL$ENSEMBL, ":", 2)[,1]

# Mapear a ENTREZID
entrez_map <- clusterProfiler::bitr(
  de_res_HL$ENSEMBL,
  fromType = "ENSEMBL",
  toType = c("ENTREZID","SYMBOL"),
  OrgDb = org.Hs.eg.db
)

# Crear nuevas columnas en de_res_HL con symbol y entrez
de_res_HL$ENTREZID <- entrez_map$ENTREZID[match(de_res_HL$ENSEMBL, entrez_map$ENSEMBL)]
de_res_HL$SYMBOL <- entrez_map$SYMBOL[match(de_res_HL$ENSEMBL, entrez_map$ENSEMBL)]

write.csv(de_res_HL, file = "norm_low_vs_high", row.names = TRUE)

#TABLA COMP IL : extraemos conteos x muestra
count_cols <- grep("^countData", colnames(IL_df))
count_mat <- IL_df[, count_cols]
rownames(count_mat) <- IL_df$X

#empezamos el analisis con la nueva matriz
DGEList<- DGEList(count_mat, remove.zeros = TRUE)
groups <- c("interm", "interm", "low", "interm", "low", "low")
clinical_status <- data.frame(
  sample = colnames(count_mat), 
  group = factor(groups, levels = c("low", "interm", "high"))
)
design <- model.matrix(~0 + group, data = clinical_status)

keep.expr<-filterByExpr(DGEList,design)
dge<-DGEList[keep.expr,]

log_counts <- log2(dge$counts + 1)

dge<- calcNormFactors(dge, method = "TMM")
tmmexp<-cpm(dge, log = T, prior.count = 3)

fit<-lmFit(tmmexp,design)

cont_mat <- makeContrasts(
  high_vs_low = groupinterm - grouplow,
  
  levels = design
)
fit2_IL <- contrasts.fit(fit, cont_mat)
fit2_IL <- eBayes(fit2_IL)

de_res_IL <- topTable(fit2_IL, number = Inf)
de_res_IL$ENSEMBL <- rownames(de_res_IL)

de_res_IL$ENSEMBL<- stringr::str_split_fixed(de_res_IL$ENSEMBL, ":", 2)[,1]

# Mapear a ENTREZID
entrez_map <- clusterProfiler::bitr(
  de_res_IL$ENSEMBL,
  fromType = "ENSEMBL",
  toType = c("ENTREZID","SYMBOL"),
  OrgDb = org.Hs.eg.db
)

# Crear nuevas columnas en de_res_IL con symbol y entrez
de_res_IL$ENTREZID <- entrez_map$ENTREZID[match(de_res_IL$ENSEMBL, entrez_map$ENSEMBL)]
de_res_IL$SYMBOL <- entrez_map$SYMBOL[match(de_res_IL$ENSEMBL, entrez_map$ENSEMBL)]

write.csv(de_res_IL, file = "norm_iterm_vs_low", row.names = TRUE)

#TABLA COMP HI : extraemos conteos x muestra
count_cols <- grep("^countData", colnames(HI_df))
count_mat <- HI_df[, count_cols]
rownames(count_mat) <- HI_df$X

#empezamos el analisis con la nueva matriz
DGEList<- DGEList(count_mat, remove.zeros = TRUE)
groups <- c("interm", "high", "interm", "interm", "high", "high")
clinical_status <- data.frame(
  sample = colnames(count_mat), 
  group = factor(groups, levels = c("low", "interm", "high"))
)
design <- model.matrix(~0 + group, data = clinical_status)

keep.expr<-filterByExpr(DGEList,design)
dge<-DGEList[keep.expr,]

log_counts <- log2(dge$counts + 1)

dge<- calcNormFactors(dge, method = "TMM")
tmmexp<-cpm(dge, log = T, prior.count = 3)
fit<-lmFit(tmmexp,design)

cont_mat <- makeContrasts(
  high_vs_low = grouphigh - groupinterm,
  
  levels = design
)

fit2_HI <- contrasts.fit(fit, cont_mat)
fit2_HI <- eBayes(fit2_HI)

de_res_HI <- topTable(fit2_HI, number = Inf)
de_res_HI$ENSEMBL <- rownames(de_res_HI)
de_res_HI$ENSEMBL<- stringr::str_split_fixed(de_res_HI$ENSEMBL, ":", 2)[,1]

# Mapear a ENTREZID
entrez_map <- clusterProfiler::bitr(
  de_res_HI$ENSEMBL,
  fromType = "ENSEMBL",
  toType = c("ENTREZID","SYMBOL"),
  OrgDb = org.Hs.eg.db
)

# Crear nuevas columnas en de_res_HI con symbol y entrez
de_res_HI$ENTREZID <- entrez_map$ENTREZID[match(de_res_HI$ENSEMBL, entrez_map$ENSEMBL)]
de_res_HI$SYMBOL <- entrez_map$SYMBOL[match(de_res_HI$ENSEMBL, entrez_map$ENSEMBL)]



write.csv(de_res_HI, file = "norm_high_vs_interm", row.names = TRUE)


#exones con Padj significantes por grupos quitando el NA:
LH_sigexons <- de_res_HL[!is.na(de_res_HL$adj.P.Val) & (de_res_HL$adj.P.Val < 0.4), ]
IL_sigexons <- de_res_IL[!is.na(de_res_IL$adj.P.Val) & (de_res_IL$adj.P.Val < 0.4), ]
HI_sigexons <- de_res_HI[!is.na(de_res_HI$adj.P.Val) & (de_res_HI$adj.P.Val < 0.4), ]


#Dataframe con todos los exones significantes de todas las comparativas
all_sigexons <- bind_rows(
  LH_sigexons %>% mutate(comparative = "LH"),
  HI_sigexons %>% mutate(comparative = "HI"),
  IL_sigexons %>% mutate(comparative = "IL"))

# Si el identificador está en rownames, primero lo pasamos a columna
all_sigexons <- tibble::rownames_to_column(all_sigexons, var = "exons")

df_upset <- all_sigexons %>%
  dplyr::select(exons, comparative) %>%
  distinct() %>%              # evita duplicados
  mutate(value = 1) %>%
  tidyr::pivot_wider(
    names_from = comparative, 
    values_from = value, 
    values_fill = 0
  )
#sacamos un upset diagram
upset(
  df_upset,
  intersect = c("LH", "HI", "IL"),
  name = "Exons",
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 3)
    ) + 
      labs(y = "Conteo")
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

signif_genes_HL <- na.omit(unique(de_res_HL$SYMBOL[de_res_HL$adj.P.Val < 0.4]))
signif_genes_HI <- na.omit(unique(de_res_HI$SYMBOL[de_res_HI$adj.P.Val < 0.4]))
signif_genes_IL <- na.omit(unique(de_res_IL$SYMBOL[de_res_IL$adj.P.Val < 0.4]))
sig <- list(
  HL = signif_genes_HL,
  HI = signif_genes_HI,
  IL = signif_genes_IL
)
mat <- list_to_matrix(sig)
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




#VennDiagram con los tres grupos. Sólo exones significativos.
venn_matrix <- as.matrix(df[, c("HL", "HI", "IL")])
vennDiagram(
  venn_matrix,include= c("both") ,
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

#volcano para LH 

de_res_HL$diffexpressed <- "NO"
de_res_HL$diffexpressed[de_res_HL$logFC > 2 & de_res_HL$adj.P.Val < 0.4] <- "UP"
de_res_HL$diffexpressed[de_res_HL$logFC < -2 & de_res_HL$adj.P.Val < 0.4] <- "DOWN"

de_res_HL$delabel <- NA

mycolors <- c("DOWN" = "blue", "UP" = "red", "NO" = "black")

de_res_HL <- de_res_HL %>%
  mutate(delabel = ifelse(
    row_number() %in% order(adj.P.Val)[1:10], 
    SYMBOL,  
    NA
  ))

# Volcano plot

ggplot(de_res_HL, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed, label = delabel)) +
  geom_point() +
  scale_color_manual(values = mycolors, name = "Regulation\n(padj < 0.4 & |logFC| > 2)") +
  geom_text_repel(size = 3, max.overlaps = 10) +
  theme_minimal() +  
  labs(title = "High vs Low")

#volcano para HI 

de_res_HI$diffexpressed <- "NO"
de_res_HI$diffexpressed[de_res_HI$logFC > 2 & de_res_HI$adj.P.Val < 0.4] <- "UP"
de_res_HI$diffexpressed[de_res_HI$logFC < -2 & de_res_HI$adj.P.Val < 0.4] <- "DOWN"

de_res_HI$delabel <- NA

mycolors <- c("DOWN" = "blue", "UP" = "red", "NO" = "black")

de_res_HI <- de_res_HI %>%
  mutate(delabel = ifelse(
    row_number() %in% order(adj.P.Val)[1:10], 
    SYMBOL,  
    NA
  ))
ggplot(de_res_HI, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed, label = delabel)) +
  geom_point() +
  scale_color_manual(values = mycolors , name = "Regulation\n(padj < 0.4 & |logFC| > 2)")  + 
  geom_text_repel(size = 3, max.overlaps = 15) +
  theme_minimal() +  
  labs(title = "High vs Interm")

#volcano para IL de todos

de_res_IL$diffexpressed <- "NO"
de_res_IL$diffexpressed[de_res_IL$logFC > 2 & de_res_IL$adj.P.Val < 0.4] <- "UP"
de_res_IL$diffexpressed[de_res_IL$logFC < -2 & de_res_IL$adj.P.Val < 0.4] <- "DOWN"

de_res_IL$delabel <- NA

mycolors <- c("DOWN" = "blue", "UP" = "red", "NO" = "black")

de_res_IL <- de_res_IL %>%
  mutate(delabel = ifelse(
    row_number() %in% order(adj.P.Val)[1:10], 
    SYMBOL,  
    NA
  ))
ggplot(de_res_IL, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed, label = delabel)) +
  geom_point() +
  scale_color_manual(values = mycolors, name = "Regulation\n(padj < 0.4 & |logFC| > 2)") +
  geom_text_repel(size = 3, max.overlaps = 10) +
  theme_minimal() +  
  labs(title = "Interm vs Low")

#HEATMAP para HL

