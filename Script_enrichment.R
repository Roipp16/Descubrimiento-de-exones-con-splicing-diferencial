library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(clusterProfiler)
library(ggplot2)
library(DOSE)


#HIGH INTERM ANALYSIS
sig_genes_HI<- read.csv("./norm_high_vs_interm")
go_genes_hi<- sig_genes_HI[sig_genes_HI$P.Value< 0.05, ]
gene_ids_hi <- go_genes_hi$ENTREZID[!is.na(go_genes_hi$ENTREZID)]

ego_hi <- enrichGO(gene         = gene_ids_hi,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",         
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

#goplot(ego)    
df_ego_hi <- ego_hi@result

top_terms_hi <- df_ego_hi[1:min(20, nrow(df_ego_hi)), ]

ggplot(top_terms_hi, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "GO Enrichment High vs Interm", x = "GO Term", y = "Gene Count") +
  theme_minimal()

KEGG_genes_hi <- enrichKEGG(gene         = gene_ids_hi, organism     = "hsa", pvalueCutoff = 0.05)
df_kegg_hi <- KEGG_genes_hi@result
top_kegg_hi <- df_kegg_hi[1:min(20, nrow(df_kegg_hi)), ]

ggplot(top_kegg_hi, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "KEGG High vs Interm",
       x = "KEGG Pathway",
       y = "Gene Count",
       fill = "p-value") +
  theme_minimal()


gene_stat_named_hi <- go_genes_hi %>%
  filter(!is.na(ENTREZID)) %>%            
  group_by(ENTREZID) %>%
  summarise(logFC = max(logFC), .groups="drop")

gene_stat_named_hi <- setNames(gene_stat_named_hi$logFC, gene_stat_named_hi$ENTREZID)
gene_stat_named_hi <- sort(gene_stat_named_hi, decreasing = TRUE)

ego2 <- gseGO(
  geneList      = gene_stat_named_hi,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  minGSSize     = 100,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

dotplot(ego2, 
        showCategory = 10, 
        title = "GSE High vs Interm",
        color = "pvalue")

ridgeplot(ego2) + 
  ggtitle("Ridge High vs Interm") +  # Add title
  theme(
    axis.text.y = element_text(family = "Courier", size = 12),  # Customize y-axis text
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)  # Center & style title
  )

#INTERM VS LOW ANALYSIS

sig_genes_IL<- read.csv("./norm_iterm_vs_low")
go_genes_IL<- sig_genes_IL[sig_genes_IL$P.Value< 0.05, ]
gene_ids_IL <- go_genes_IL$ENTREZID[!is.na(go_genes_IL$ENTREZID)]

ego_IL <- enrichGO(gene         = gene_ids_IL,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "BP",         
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

#goplot(ego)    
df_ego_IL <- ego_IL@result

top_terms_IL <- df_ego_IL[1:min(20, nrow(df_ego_IL)), ]

ggplot(top_terms_IL, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "GO Enrichment Interm vs Low", x = "GO Term", y = "Gene Count") +
  theme_minimal()

KEGG_genes_IL <- enrichKEGG(gene         = gene_ids_IL, organism     = "hsa", pvalueCutoff = 0.05)
df_kegg_IL <- KEGG_genes_IL@result
top_kegg_IL <- df_kegg_IL[1:min(20, nrow(df_kegg_IL)), ]

ggplot(top_kegg_IL, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "KEGG Interm vs Low",
       x = "KEGG Pathway",
       y = "Gene Count",
       fill = "p-value") +
  theme_minimal()


gene_stat_named_IL <- go_genes_IL %>%
  filter(!is.na(ENTREZID)) %>%            
  group_by(ENTREZID) %>%
  summarise(logFC = max(logFC), .groups="drop")

gene_stat_named_IL <- setNames(gene_stat_named_IL$logFC, gene_stat_named_IL$ENTREZID)
gene_stat_named_IL <- sort(gene_stat_named_IL, decreasing = TRUE)

ego2 <- gseGO(
  geneList      = gene_stat_named_IL,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  minGSSize     = 100,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

dotplot(ego2, 
        showCategory = 10, 
        title = "GSE Interm vs Low",
        color = "pvalue")

ridgeplot(ego2) + 
  ggtitle("Ridge Interm vs Low") + 
  theme(
    axis.text.y = element_text(family = "Courier", size = 12),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14) 
  )



#HIGH VS LOW ANALYSIS


sig_genes_HL<- read.csv("./norm_low_vs_high")
go_genes_HL<- sig_genes_HL[sig_genes_HL$P.Value< 0.05, ]
gene_ids_HL <- go_genes_HL$ENTREZID[!is.na(go_genes_HL$ENTREZID)]

ego_HL <- enrichGO(gene         = gene_ids_HL,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "BP",         
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

#goplot(ego)    
df_ego_HL <- ego_HL@result

top_terms_HL <- df_ego_HL[1:min(20, nrow(df_ego_HL)), ]

ggplot(top_terms_HL, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "GO Enrichment Low vs High", x = "GO Term", y = "Gene Count") +
  theme_minimal()

KEGG_genes_HL <- enrichKEGG(gene         = gene_ids_HL, organism     = "hsa", pvalueCutoff = 0.05)
df_kegg_HL <- KEGG_genes_HL@result
top_kegg_HL <- df_kegg_HL[1:min(20, nrow(df_kegg_HL)), ]

ggplot(top_kegg_HL, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "KEGG Low vs High",
       x = "KEGG Pathway",
       y = "Gene Count",
       fill = "p-value") +
  theme_minimal()


gene_stat_named_HL <- go_genes_HL %>%
  filter(!is.na(ENTREZID)) %>%            
  group_by(ENTREZID) %>%
  summarise(logFC = max(logFC), .groups="drop")

gene_stat_named_HL <- setNames(gene_stat_named_HL$logFC, gene_stat_named_HL$ENTREZID)
gene_stat_named_HL <- sort(gene_stat_named_HL, decreasing = TRUE)

ego2 <- gseGO(
  geneList      = gene_stat_named_HL,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  minGSSize     = 100,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

dotplot(ego2, 
        showCategory = 10, 
        title = "GSE Low vs High",
        color = "pvalue")

ridgeplot(ego2) + 
  ggtitle("Ridge Low vs High") + 
  theme(
    axis.text.y = element_text(family = "Courier", size = 12),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14) 
  )
