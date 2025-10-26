library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(clusterProfiler)
library(ggplot2)
library(DOSE)

#HIGH INTERM ANALYSIS
sig_genes_HI<- read.csv("~/Documents/Analysis/DEXSEQ_high_vs_interm.csv")
go_genes<- sig_genes_HI[sig_genes_HI$P.Value < 0.05, ]
gene_symbols <- go_genes$X

#haremos para los genes especificados tambien: 
selected_HI <- genes[genes %in% gene_symbols] #aqui los tres estan presentes


entrez_genes <- bitr(gene_symbols,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene         = entrez_genes$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",         
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)
          
#goplot(ego)    
df_ego <- ego@result

top_terms <- df_ego[1:min(20, nrow(df_ego)), ]

ggplot(top_terms, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "GO Enrichment High vs Interm", x = "GO Term", y = "Gene Count") +
  theme_minimal()
  
KEGG_genes <- enrichKEGG(gene         = entrez_genes$ENTREZID, organism     = "hsa", pvalueCutoff = 0.05)
df_kegg <- KEGG_genes@result
top_kegg <- df_kegg[1:min(20, nrow(df_kegg)), ]

ggplot(top_kegg, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "KEGG HI",
       x = "KEGG Pathway",
       y = "Gene Count",
       fill = "p-value") +
  theme_minimal()


gene_stat <- go_genes$logFC
names(gene_stat) <- go_genes$X  
entrez_map <- bitr(names(gene_stat), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_stat_named <- gene_stat[entrez_map$SYMBOL]
names(gene_stat_named) <- entrez_map$ENTREZID
gene_stat_named <- sort(gene_stat_named, decreasing = TRUE)

ego2 <- gseGO(geneList      = gene_stat_named,
              OrgDb         = org.Hs.eg.db,
              ont           = "BP",
              minGSSize     = 100,
              maxGSSize     = 500,
              pvalueCutoff  = 0.05,
              verbose       = FALSE)

dotplot(ego2, 
        showCategory = 10, 
        title = "GSEA HI",
        color = "pvalue")

ridgeplot(ego2) + 
  ggtitle("Ridge HI") +  # Add title
  theme(
    axis.text.y = element_text(family = "Courier", size = 12),  # Customize y-axis text
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)  # Center & style title
  )

#INTERM VS LOW ANALYSIS

sig_genes_IL<- read.csv("~/Documents/Analysis/DEXSEQ_interm_vs_low.csv")
go_genes<- sig_genes_IL[sig_genes_IL$P.Value < 0.05, ]
gene_symbols <- go_genes$X
selected_IL <- genes[genes %in% gene_symbols] #"B4GALT6" "ST8SIA3", solo estos dos

entrez_genes <- bitr(gene_symbols,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene         = entrez_genes$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",         
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

#goplot(ego)    
df_ego <- ego@result

top_terms <- df_ego[1:min(20, nrow(df_ego)), ]

ggplot(top_terms, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "GO IL", x = "GO Term", y = "Gene Count") +
  theme_minimal()

KEGG_genes <- enrichKEGG(gene         = entrez_genes$ENTREZID, organism     = "hsa", pvalueCutoff = 0.05)
df_kegg <- KEGG_genes@result
top_kegg <- df_kegg[1:min(20, nrow(df_kegg)), ]

ggplot(top_kegg, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "KEGG IL",
       x = "KEGG Pathway",
       y = "Gene Count",
       fill = "p-value") +
  theme_minimal()


gene_stat <- go_genes$logFC
names(gene_stat) <- go_genes$X  
entrez_map <- bitr(names(gene_stat), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_stat_named <- gene_stat[entrez_map$SYMBOL]
names(gene_stat_named) <- entrez_map$ENTREZID
gene_stat_named <- sort(gene_stat_named, decreasing = TRUE)

ego2 <- gseGO(geneList      = gene_stat_named,
              OrgDb         = org.Hs.eg.db,
              ont           = "BP",
              minGSSize     = 100,
              maxGSSize     = 500,
              pvalueCutoff  = 0.05,
              verbose       = FALSE)

dotplot(ego2, 
        showCategory = 10, 
        title = "GSEA IL",
        color = "pvalue")

ridgeplot(ego2) + 
  ggtitle("Ridge IL") +  # Add title
  theme(
    axis.text.y = element_text(family = "Courier", size = 12),  # Customize y-axis text
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)  # Center & style title
  )


#HIGH VS LOW ANALYSIS

sig_genes_HL<- read.csv("~/Documents/Analysis/DEXSEQ_low_vs_high.csv")
go_genes<- sig_genes_HL[sig_genes_HL$P.Value < 0.05, ]
gene_symbols <- go_genes$X
selected_HL <- genes[genes %in% gene_symbols] #ninguno

entrez_genes <- bitr(gene_symbols,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene         = entrez_genes$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",         
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

#goplot(ego)    
df_ego <- ego@result

top_terms <- df_ego[1:min(20, nrow(df_ego)), ]

ggplot(top_terms, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "GO HL", x = "GO Term", y = "Gene Count") +
  theme_minimal()

KEGG_genes <- enrichKEGG(gene         = entrez_genes$ENTREZID, organism     = "hsa", pvalueCutoff = 0.05)
df_kegg <- KEGG_genes@result
top_kegg <- df_kegg[1:min(20, nrow(df_kegg)), ]

ggplot(top_kegg, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "KEGG HL",
       x = "KEGG Pathway",
       y = "Gene Count",
       fill = "p-value") +
  theme_minimal()


gene_stat <- go_genes$logFC
names(gene_stat) <- go_genes$X  
entrez_map <- bitr(names(gene_stat), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_stat_named <- gene_stat[entrez_map$SYMBOL]
names(gene_stat_named) <- entrez_map$ENTREZID
gene_stat_named <- sort(gene_stat_named, decreasing = TRUE)

ego2 <- gseGO(geneList      = gene_stat_named,
              OrgDb         = org.Hs.eg.db,
              ont           = "BP",
              minGSSize     = 100,
              maxGSSize     = 500,
              pvalueCutoff  = 0.05,
              verbose       = FALSE)

dotplot(ego2, 
        showCategory = 10, 
        title = "GSEA HL",
        color = "pvalue")

ridgeplot(ego2) + 
  ggtitle("Ridge HL") +  # Add title
  theme(
    axis.text.y = element_text(family = "Courier", size = 12),  # Customize y-axis text
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)  # Center & style title
  )
