
#' Microarray analysis 
# GSE76039-- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76039
# Paper : https://www.jci.org/articles/view/85271#SEC4 
#package required for microarray !! 
BiocManager::install(c("oligo", "limma", "Biobase", "pd.hg.u133.plus.2"))  

#Analysis workflow for CEL files in R

#Install required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("oligo", "limma", "Biobase", "pd.hg.u133.plus.2"))  
# Install required packages
BiocManager::install(c("affy", "limma", "hgu133plus2.db", "annotate"))

#read tar file 
untar("GSE76039_RAW.tar", exdir = "GSE76039_tar_R")

# Load libraries
library(affy)          # for reading CEL files (3’ expression arrays)
library(limma)         # for differential expression
library(hgu133plus2.db) # annotation package for HG-U133 Plus 2.0
library(annotate)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)

# Read CEL files
# Set path to CEL files
cel_path <- "./ATC_Microarray/"   

# List all CEL files
celFiles <- list.celfiles(cel_path, full.names = TRUE)

# Read CEL files
rawData <- ReadAffy(filenames = celFiles)
view(rawData$sample)

# 3. Normalize (RMA)
# Robust Multi-array Average (RMA) normalization
normData <- rma(rawData)

# Extract expression matrix
exprs_matrix <- exprs(normData)
dim(exprs_matrix)   # probes x samples

# 4. Add gene symbols
# Map probe IDs to gene symbols
probe_ids <- rownames(exprs_matrix)
gene_symbols <- getSYMBOL(probe_ids, "hgu133plus2.db")

# Remove probes without gene symbols
keep <- !is.na(gene_symbols)
exprs_matrix <- exprs_matrix[keep, ]
gene_symbols <- gene_symbols[keep]

# Add annotation
exprs_annot <- data.frame(GeneSymbol = gene_symbols, exprs_matrix)

# Collapse to one row per gene (average across probes)
exprs_gene <- avereps(exprs_matrix, ID = gene_symbols)
dim(exprs_gene)  # genes x samples
view(exprs_gene)

#write csv file 
ATC_microarray_normalizedcount<- read.csv("Normalized.count_ATC.csv")

#round it to integer
df_rounded <- ATC_microarray_normalizedcount %>%
  mutate(across(where(is.numeric), round, 0))

# View the result
head(df_rounded)

# Optional: save to file
write.csv(df_rounded, "ATC_Normalized_count.csv", row.names = FALSE)
write.csv(exprs_gene, file= "Normalized.count_ATC.csv")

#load metadata
metadata <- read.csv("metadata.csv", 
                     stringsAsFactors = TRUE)

# Ensure order matches columns of exprs_matrix
metadata <- metadata[match(colnames(exprs_matrix), metadata$Microarray_long_id), ]
all(colnames(exprs_matrix) == metadata$Microarray_long_id)


#Principal component analysis
# Transpose so samples are rows
exprs_t <- t(exprs_gene)  # now samples x genes

# Run PCA (scale = TRUE standardizes each gene)
pca_res <- prcomp(exprs_t, scale. = TRUE)

# Prepare a data frame for plotting
pca_df <- data.frame(pca_res$x)
pca_df$Group <- metadata$Group  # replace 'Group' with your column of interest

# Check variance explained
summary(pca_res)

# Plot PCA
library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  theme(plot.title = element_text(hjust= 0.5)) +
  labs(
    title = "Principal component analysis",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")
  )

#' Differential gene expression =======================================================================
design <- model.matrix(~ 0 + factor(metadata$Group))
colnames(design) <- levels(factor(metadata$Group))

# Contrast male vs female
contrast.matrix <- makeContrasts(Young_M_vs_Young_F = Young_M - Young_F, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)
write.csv(results, file= "DEG_Young_male.vs_young_female.csv")


# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Convert Significance to factor to fix legend issues
all_genes$Significance <- factor(all_genes$Significance, 
                                 levels = c("Upregulated", "Not significant", "Downregulated"))

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
Young_M_vs_Young_F <- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "#BB0C00",
                                "Not significant" = "grey",
                                "Downregulated" = "#00AFBB")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Young_Male vs Young_female") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, 
                  aes(x = logFC, y = logP, label = rownames(top_genes)),
                  size = 5, 
                  box.padding = 0.3, 
                  point.padding = 0.2,
                  max.overlaps = 20)
Young_M_vs_Young_F 

# Save with high resolution
ggsave(filename = "Young_Male_vs_Young_female.tiff",    
       plot = Young_M_vs_Young_F,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in")  

#Classic enrichment 
gene_ranks <- results$t
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)














# download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/ 
gmt_file <- "h.all.v7.5.1.symbols.gmt"
gene_sets <- read.gmt(gmt_file)

# Run GSEA
gsea_results <- GSEA(
  gene_ranks,
  TERM2GENE = gene_sets,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View results
head(gsea_results@result)
gsea_table <- as.data.frame(gsea_results)
write.csv(gsea_table, file= "Hallmark_Young_male_vs_young_female.csv")

# Add sign column==== very useful code 
hallmark_df <- as.data.frame(gsea_results)
hallmark_df$sign <- ifelse(hallmark_df$NES > 0, "Young male", "Young female")
gsea_results@result$sign <- hallmark_df$sign

# Dotplot split by NES sign
hallmark.dotplot<- dotplot(gsea_results, showCategory = 15, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hallmark Pathways")
hallmark.dotplot

#hallmark
ggsave(filename = "hallmark_Young_male_vs_Young_female_all.tiff",    
       plot = hallmark.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


#KEGG Pathway
install.packages("msigdbr") # if not installed
library(msigdbr)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)

# Use t-statistics or logFC as ranking metric
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

# Map gene symbols to Entrez IDs
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# Keep only mapped genes
gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID

# Sort again just in case
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

## Run KEGG GSEA 
kegg_gsea <- gseKEGG(
  geneList     = gene_ranks,
  organism     = "hsa",          # human
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# View results
head(kegg_gsea@result)
kegg_table <- as.data.frame(kegg_gsea)
write.csv(kegg_table, file= "KEGG_Young male_vs_young female.csv")

# Add sign column==== very useful code 
kegg_df <- as.data.frame(kegg_gsea)
kegg_df$sign <- ifelse(kegg_df$NES > 0, "Young male", "Young female")
kegg_gsea@result$sign <- kegg_df$sign

# Dotplot split by NES sign
KEGG.dotplot<- dotplot(kegg_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("KEGG Pathways")
KEGG.dotplot

#save plot
ggsave(filename = "KEGG_young_male_vs_Young_female_all.tiff",    
       plot = KEGG.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


# Gene Ontology 
library(clusterProfiler)
library(org.Hs.eg.db)

# ----- Gene Ranks -----
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Map SYMBOL -> ENTREZ
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

#  Run GO GSEA 
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_BP_Young_male_vs_Young_female.csv")

# Add sign column (male = NES > 0, female = NES < 0) 
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "Young male", "Young female")
go_gsea@result$sign <- go_df$sign

# ----- Dotplot split by NES sign -----
GO.dotplot <- dotplot(go_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Biological Process Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_Young_male_vs_Young_female_all.tiff",    
       plot = GO.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


# GO molecular function 
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_Molecular function_Young_male_vs_Young_female.csv")

# Add sign column (male = NES > 0, female = NES < 0) 
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "Young male", "Young female")
go_gsea@result$sign <- go_df$sign

# ----- Dotplot split by NES sign -----
GO.dotplot <- dotplot(go_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Biological Process Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_Young_male_vs_Young_female_all.tiff",    
       plot = GO.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


#=============================================================================================
                           ## old male vs. old female
#============================================================================================
# Differential gene expression 
design <- model.matrix(~ 0 + factor(metadata$Group))
colnames(design) <- levels(factor(metadata$Group))

# Contrast male vs female
contrast.matrix <- makeContrasts(Old_M_vs_Old_F = Old_M - Old_F, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)
write.csv(results, file= "DEG_Old_male.vs_Old_female.csv")

library(ggplot2)
library(ggrepel)

# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Convert Significance to factor to fix legend issues
all_genes$Significance <- factor(all_genes$Significance, 
                                 levels = c("Upregulated", "Not significant", "Downregulated"))

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
Old_M_vs_Old_F <- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "#BB0C00",
                                "Not significant" = "grey",
                                "Downregulated" = "#00AFBB")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Old_Male vs Old_female") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, 
                  aes(x = logFC, y = logP, label = rownames(top_genes)),
                  size = 5, 
                  box.padding = 0.3, 
                  point.padding = 0.2,
                  max.overlaps = 20)
Old_M_vs_Old_F 

# Save with high resolution
ggsave(filename = "Old_Male_vs_Old_female.tiff",    
       plot = Old_M_vs_Old_F ,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in")  

#Classic enrichment 
gene_ranks <- results$t
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

# run gene set enrichment analalysis via cluster profiler
# download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/ 
gmt_file <- "h.all.v7.5.1.symbols.gmt"
gene_sets <- read.gmt(gmt_file)

# Run GSEA
gsea_results <- GSEA(
  gene_ranks,
  TERM2GENE = gene_sets,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View results
head(gsea_results@result)
gsea_table <- as.data.frame(gsea_results)
write.csv(gsea_table, file= "Hallmark_Old_male_vs_Old_female.csv")

# Add sign column==== very useful code 
hallmark_df <- as.data.frame(gsea_results)
hallmark_df$sign <- ifelse(hallmark_df$NES > 0, "Old male", "Old female")
gsea_results@result$sign <- hallmark_df$sign

# Dotplot split by NES sign
hallmark.dotplot<- dotplot(gsea_results, showCategory = 15, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hallmark Pathways")
hallmark.dotplot

#hallmark
ggsave(filename = "hallmark_old_male_vs_old_female_all.tiff",    
       plot = hallmark.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")

#' KEGG Pathway ===========================================================================================
install.packages("msigdbr") # if not installed
library(msigdbr)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)

# Use t-statistics or logFC as ranking metric
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

# Map gene symbols to Entrez IDs
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# Keep only mapped genes
gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID

# Sort again just in case
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

## Run KEGG GSEA 
kegg_gsea <- gseKEGG(
  geneList     = gene_ranks,
  organism     = "hsa",          # human
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# View results
head(kegg_gsea@result)
kegg_table <- as.data.frame(kegg_gsea)
write.csv(kegg_table, file= "KEGG_Old male_vs_Old female.csv")

# Add sign column==== very useful code 
kegg_df <- as.data.frame(kegg_gsea)
kegg_df$sign <- ifelse(kegg_df$NES > 0, "Old male", "Old female")
kegg_gsea@result$sign <- kegg_df$sign

# Dotplot split by NES sign
KEGG.dotplot<- dotplot(kegg_gsea, showCategory = 20, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("KEGG Pathways")
KEGG.dotplot

#save plot
ggsave(filename = "KEGG_old_male_vs_old_female_all.tiff",    
       plot = KEGG.dotplot,                        
       width = 7, 
       height = 9,           
       dpi = 600,                          
       units = "in")


# Gene Ontology 
library(clusterProfiler)
library(org.Hs.eg.db)

#  Gene Ranks 
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Map SYMBOL -> ENTREZ
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

#  Run GO GSEA 
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_BP_Old_male_vs_Old_female.csv")

# Add sign column (male = NES > 0, female = NES < 0) 
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "Old male", "Old female")
go_gsea@result$sign <- go_df$sign

#' Dotplot split by NES sign 
GO.dotplot <- dotplot(go_gsea, showCategory = 15, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Biological Process Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_BP_old_male_vs_old_female_all.tiff",    
       plot = GO.dotplot,                        
       width = 7, 
       height = 9,           
       dpi = 600,                          
       units = "in")


# GO molecular function 
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_Molecular function_Young_male_vs_Young_female.csv")

# Add sign column (male = NES > 0, female = NES < 0) 
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "Young male", "Young female")
go_gsea@result$sign <- go_df$sign

#  Dotplot split by NES sign 
GO.dotplot <- dotplot(go_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Biological Process Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_Young_male_vs_Young_female_all.tiff",    
       plot = GO.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")

#======================================================================================
                                ## Young famele vs. Old female
#======================================================================================

# Differential gene expression 
design <- model.matrix(~ 0 + factor(metadata$Group))
colnames(design) <- levels(factor(metadata$Group))

# Contrast male vs female
contrast.matrix <- makeContrasts(Old_F_vs_young_F = Old_F - Young_F, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)
write.csv(results, file= "DEG_Old_female.vs_young_female.csv")

library(ggplot2)
library(ggrepel)

# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Convert Significance to factor to fix legend issues
all_genes$Significance <- factor(all_genes$Significance, 
                                 levels = c("Upregulated", "Not significant", "Downregulated"))

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
Old_F_vs_young_F <- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "#BB0C00",
                                "Not significant" = "grey",
                                "Downregulated" = "#00AFBB")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Old_Female vs Young_female") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, 
                  aes(x = logFC, y = logP, label = rownames(top_genes)),
                  size = 5, 
                  box.padding = 0.3, 
                  point.padding = 0.2,
                  max.overlaps = 20)
Old_F_vs_young_F

# Save with high resolution
ggsave(filename = "Old_Female_vs_Young_female.tiff",    
       plot = Old_F_vs_young_F ,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in")  

#Classic enrichment 
gene_ranks <- results$t
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

# download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/ 
gmt_file <- "h.all.v7.5.1.symbols.gmt"
gene_sets <- read.gmt(gmt_file)

# Run GSEA
gsea_results <- GSEA(
  gene_ranks,
  TERM2GENE = gene_sets,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View results
head(gsea_results@result)
gsea_table <- as.data.frame(gsea_results)
write.csv(gsea_table, file= "Hallmark_Old_female_vs_young_female.csv")

# Add sign column==== very useful code 
hallmark_df <- as.data.frame(gsea_results)
hallmark_df$sign <- ifelse(hallmark_df$NES > 0, "Old female", "Young female")
gsea_results@result$sign <- hallmark_df$sign

# Dotplot split by NES sign
hallmark.dotplot<- dotplot(gsea_results, showCategory = 15, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hallmark Pathways")
hallmark.dotplot

#hallmark
ggsave(filename = "hallmark_old_female_vs_young_female_all.tiff",    
       plot = hallmark.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")

#' KEGG Pathway ====================================================================================================
install.packages("msigdbr") # if not installed
library(msigdbr)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)

# Use t-statistics or logFC as ranking metric
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

# Map gene symbols to Entrez IDs
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# Keep only mapped genes
gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID

# Sort again just in case
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

## Run KEGG GSEA 
kegg_gsea <- gseKEGG(
  geneList     = gene_ranks,
  organism     = "hsa",          # human
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# View results
head(kegg_gsea@result)
kegg_table <- as.data.frame(kegg_gsea)
write.csv(kegg_table, file= "KEGG_Old female_vs_young female.csv")

# Add sign column==== very useful code 
kegg_df <- as.data.frame(kegg_gsea)
kegg_df$sign <- ifelse(kegg_df$NES > 0, "Old female", "Young female")
kegg_gsea@result$sign <- kegg_df$sign

# Dotplot split by NES sign
KEGG.dotplot<- dotplot(kegg_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("KEGG Pathways")
KEGG.dotplot


#save plot
ggsave(filename = "KEGG_old_female_vs_young_female_all.tiff",    
       plot = KEGG.dotplot,                        
       width = 7, 
       height = 9,           
       dpi = 600,                          
       units = "in")


# Gene Ontology 
library(clusterProfiler)
library(org.Hs.eg.db)

# ----- Gene Ranks -----
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Map SYMBOL -> ENTREZ
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

#  Run GO GSEA 
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_BP_Old_female_vs_young_female.csv")

# Add sign column (male = NES > 0, female = NES < 0) 
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "Old female", "Young female")
go_gsea@result$sign <- go_df$sign

#'  Dotplot split by NES sign ============================================================================
GO.dotplot <- dotplot(go_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Biological Process Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_BP_old_female_vs_young_female_all.tiff",    
       plot = GO.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


# GO molecular function 
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_Molecular function_Young_female_vs_Young_female.csv")

# Add sign column (male = NES > 0, female = NES < 0) 
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "old female", "Young female")
go_gsea@result$sign <- go_df$sign

# ----- Dotplot split by NES sign -----
GO.dotplot <- dotplot(go_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Melecular function Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_MF_Old female_vs_Young_female_all.tiff",    
       plot = GO.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")

##========================================================================================= 
                                 ## male old vs. Male young 
##=========================================================================================
# Differential gene expression 
design <- model.matrix(~ 0 + factor(metadata$Group))
colnames(design) <- levels(factor(metadata$Group))

# Contrast male vs female
contrast.matrix <- makeContrasts(Old_F_vs_young_F = Old_M - Young_M, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)
write.csv(results, file= "DEG_Old_male.vs_young_male.csv")

library(ggplot2)
library(ggrepel)

# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Convert Significance to factor to fix legend issues
all_genes$Significance <- factor(all_genes$Significance, 
                                 levels = c("Upregulated", "Not significant", "Downregulated"))

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
Old_male_vs_young_male <- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "#BB0C00",
                                "Not significant" = "grey",
                                "Downregulated" = "#00AFBB")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Old_male vs Young_male") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, 
                  aes(x = logFC, y = logP, label = rownames(top_genes)),
                  size = 5, 
                  box.padding = 0.3, 
                  point.padding = 0.2,
                  max.overlaps = 20)
Old_male_vs_young_male

# Save with high resolution
ggsave(filename = "Old_male_vs_Young_male.tiff",    
       plot = Old_male_vs_young_male ,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in")  

#Classic enrichment 
gene_ranks <- results$t
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

# download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/ 
gmt_file <- "h.all.v7.5.1.symbols.gmt"
gene_sets <- read.gmt(gmt_file)

# Run GSEA
gsea_results <- GSEA(
  gene_ranks,
  TERM2GENE = gene_sets,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View results
head(gsea_results@result)
gsea_table <- as.data.frame(gsea_results)
write.csv(gsea_table, file= "Hallmark_Old_male_vs_young_male.csv")

# Add sign column==== very useful code 
hallmark_df <- as.data.frame(gsea_results)
hallmark_df$sign <- ifelse(hallmark_df$NES > 0, "Old male", "Young male")
gsea_results@result$sign <- hallmark_df$sign

# Dotplot split by NES sign
hallmark.dotplot<- dotplot(gsea_results, showCategory = 15, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hallmark Pathways")
hallmark.dotplot

#hallmark
ggsave(filename = "hallmark_old_male_vs_young_male.tiff",    
       plot = hallmark.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")

#KEGG Pathway
install.packages("msigdbr") # if not installed
library(msigdbr)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)

# Use t-statistics or logFC as ranking metric
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

# Map gene symbols to Entrez IDs
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# Keep only mapped genes
gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID

# Sort again just in case
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

## Run KEGG GSEA 
kegg_gsea <- gseKEGG(
  geneList     = gene_ranks,
  organism     = "hsa",          # human
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# View results
head(kegg_gsea@result)
kegg_table <- as.data.frame(kegg_gsea)
write.csv(kegg_table, file= "KEGG_Old male_vs_young male.csv")

# Add sign column==== very useful code 
kegg_df <- as.data.frame(kegg_gsea)
kegg_df$sign <- ifelse(kegg_df$NES > 0, "Old male", "Young male")
kegg_gsea@result$sign <- kegg_df$sign

# Dotplot split by NES sign
KEGG.dotplot<- dotplot(kegg_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("KEGG Pathways")
KEGG.dotplot

#save plot
ggsave(filename = "KEGG_old_male_vs_young_male.tiff",    
       plot = KEGG.dotplot,                        
       width = 7, 
       height = 8,           
       dpi = 600,                          
       units = "in")

# Gene Ontology 
library(clusterProfiler)
library(org.Hs.eg.db)

#  Gene Ranks 
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Map SYMBOL -> ENTREZ
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

#  Run GO GSEA 
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_BP_Old_female_vs_young_female.csv")

# Add sign column (male = NES > 0, female = NES < 0) 
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "Old male", "Young male")
go_gsea@result$sign <- go_df$sign

#  Dotplot split by NES sign 
GO.dotplot <- dotplot(go_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Biological Process Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_BP_old_male_vs_young_male.tiff",    
       plot = GO.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


# GO molecular function 
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_Molecular function_old male_vs_Young_male.csv")

# Add sign column (male = NES > 0, female = NES < 0) 
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "old male", "Young female")
go_gsea@result$sign <- go_df$sign

# ----- Dotplot split by NES sign -----
GO.dotplot <- dotplot(go_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Melecular function Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_MF_Old male_vs_Young_male.tiff",    
       plot = GO.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


#. Differential Expression using limma
# Design matrix
design <- model.matrix(~ 0 + factor(metadata$Gender))
colnames(design) <- levels(factor(metadata$Gender))

# Contrast male vs female
contrast.matrix <- makeContrasts(male_vs_female = male - female, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)

# Filter by adjusted p-value < 0.05
deg_filtered <- subset(results, adj.P.Val < 0.05)

# Check top rows
head(deg_filtered)

# Save to CSV
write.csv(deg_filtered, "DEG_male_vs_female_FDR0.05.csv", row.names = TRUE)

## Volcano plot 
library(ggplot2)
library(ggrepel)

# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Convert Significance to factor to fix legend issues
all_genes$Significance <- factor(all_genes$Significance, 
                                 levels = c("Upregulated", "Not significant", "Downregulated"))

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
Male_vs_female_all <- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "#BB0C00",
                                "Not significant" = "grey",
                                "Downregulated" = "#00AFBB")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Male vs Female (all patients)") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, 
                  aes(x = logFC, y = logP, label = rownames(top_genes)),
                  size = 5, 
                  box.padding = 0.3, 
                  point.padding = 0.2,
                  max.overlaps = 20)
Male_vs_female_all

# Save with high resolution
ggsave(filename = "Male_vs_female.tiff",    
       plot = Male_vs_female_all,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in")  

#Classic GSEA male vs.female 
# Create ranked list for GSEA (named vector: gene → score)
BiocManager::install("clusterProfiler")
gene_ranks <- results$t
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

library(clusterProfiler)
# download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/ 
gmt_file <- "h.all.v7.5.1.symbols.gmt"
gene_sets <- read.gmt(gmt_file)

# Run GSEA
gsea_results <- GSEA(
  gene_ranks,
  TERM2GENE = gene_sets,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View results
head(gsea_results@result)
gsea_table <- as.data.frame(gsea_results)

# Add sign column==== very useful code 
hallmark_df <- as.data.frame(gsea_results)
hallmark_df$sign <- ifelse(hallmark_df$NES > 0, "male", "female")
gsea_results@result$sign <- hallmark_df$sign

# Dotplot split by NES sign
hallmark.dotplot<- dotplot(gsea_results, showCategory = 15, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hallmark Pathways")
hallmark.dotplot

#hallmark
ggsave(filename = "hallmark_male_vs_female_all.tiff",    
       plot = hallmark.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")

# ridge plot 
hallmark_male_vs_female<- ridgeplot(gsea_results) +
  theme_minimal() +
  theme(plot.title= element_text(hjust=0.5, face= "bold", size= 14),
    axis.text.y = element_text(size = 15, face= "bold"),  
    axis.text.x = element_text(size = 15, face= "bold"), 
    axis.title.x = element_text()) + 
  labs(x= "NES", title= "Male vs. Female")
hallmark_male_vs_female

#save plot
ggsave(filename = "hallmark_male_vs_female_all.tiff",    
       plot = hallmark_male_vs_female,                        
       width = 10, 
       height = 13,           
       dpi = 600,                          
       units = "in")

#visualization 
library(enrichplot)
library(ggplot2)

# Get top 5 enriched pathways (based on NES or p.adjust)
top_terms <- gsea_results@result$ID[1:5]

# Plot and save each one
for (i in seq_along(top_terms)) {
  p <- gseaplot2(gsea_results,
                 geneSetID = top_terms[i],
                 title = gsea_results@result$Description[i])
  
  # save as PDF or PNG
  ggsave(paste0("GSEA_plot_", i, ".pdf"), plot = p, width = 7, height = 5)
}

#or all result 
library(enrichplot)
# Convert GSEA results to data frame
gsea_df <- as.data.frame(gsea_results)
write.csv(gsea_df, file= "gsea_male_vs_female_allsample.csv")

# Top 5 positive NES
top_pos <- gsea_df[order(-gsea_df$NES), ][1:5, "ID"]

# Top 5 negative NES
top_neg <- gsea_df[order(gsea_df$NES), ][1:1, "ID"]

# Combine for plotting
combine <- c(top_pos, top_neg)

# Top 5 pathways
gseaplot2(gsea_results, geneSetID = 1:5, pvalue_table = TRUE)  ## use ful 

# Plot all 10 pathways in one overlay plot
top10<- gseaplot2(gsea_results,
          geneSetID = combine,
          pvalue_table = FALSE) 
top10 

ggsave(filename = "Classic_GSEA_male_vs_female.tiff",    
       plot = top10,                        
       width = 8, 
       height = 8,           
       dpi = 600,                          
       units = "in") 

#KEGG Pathway
install.packages("msigdbr") # if not installed
library(msigdbr)
library(GSEABase)

# Use t-statistics or logFC as ranking metric
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)


library(org.Hs.eg.db)
library(clusterProfiler)

# Map gene symbols to Entrez IDs
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# Keep only mapped genes
gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID

# Sort again just in case
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

## Run KEGG GSEA 
kegg_gsea <- gseKEGG(
  geneList     = gene_ranks,
  organism     = "hsa",          # human
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# Add sign column==== very useful code 
kegg_df <- as.data.frame(kegg_gsea)
kegg_df$sign <- ifelse(kegg_df$NES > 0, "male", "female")
kegg_gsea@result$sign <- kegg_df$sign

# Dotplot split by NES sign
KEGG.dotplot<- dotplot(kegg_gsea, showCategory = 10, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("KEGG Pathways")
KEGG.dotplot


#save plot
ggsave(filename = "KEGG_male_vs_female_all.tiff",    
       plot = KEGG.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


# Bar plot left and right =============
# Convert GSEA object to a data frame
gsea_df <- as.data.frame(kegg_gsea)

# Convert GSEA object to data frame
gsea_df <- as.data.frame(kegg_gsea)

# Top 10 positive NES
top_pos <- gsea_df[gsea_df$NES > 0, ]
top_pos <- top_pos[order(-top_pos$NES), ]
top_pos <- head(top_pos, 10)

# Top 10 negative NES
top_neg <- gsea_df[gsea_df$NES < 0, ]
top_neg <- top_neg[order(top_neg$NES), ]
top_neg <- head(top_neg, 10)

# Combine
top20 <- rbind(top_pos, top_neg)
top20$Description <- factor(top20$Description, levels = top20$Description[order(top20$NES)])

# Order pathways by NES for plotting  == new visualization 
ggplot(top20, aes(x = NES, y = Description, fill = NES > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#BB0C00", "FALSE" = "#00AFBB"),
                    labels = c("Negative NES", "Positive NES")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face="bold", size=12),
    axis.text.x = element_text(face="bold", size=12),
    plot.title = element_text(face="bold", size=14, hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_blank()
  ) + 
  theme(legend.position = "none") +
  labs(title = "Top 10 Positive & Top 10 Negative KEGG Pathways", x="NES", y="")

#save plot
ggsave(filename = "hallmark_male_vs_female_all.tiff",    
       plot = hallmark_male_vs_female,                        
       width = 10, 
       height = 13,           
       dpi = 600,                          
       units = "in")

#visualization 
library(enrichplot)
library(ggplot2)

# Get top 5 enriched pathways (based on NES or p.adjust)
top_terms <- gsea_results@result$ID[1:5]

# Plot and save each one
for (i in seq_along(top_terms)) {
  p <- gseaplot2(gsea_results,
                 geneSetID = top_terms[i],
                 title = gsea_results@result$Description[i])
  
  # save as PDF or PNG
  ggsave(paste0("GSEA_plot_", i, ".pdf"), plot = p, width = 7, height = 5)
}

#or all result 
library(enrichplot)
# Convert GSEA results to data frame
gsea_df <- as.data.frame(gsea_results)
write.csv(gsea_df, file= "gsea_male_vs_female_allsample.csv")

# Top 5 positive NES
top_pos <- gsea_df[order(-gsea_df$NES), ][1:5, "ID"]

# Top 5 negative NES
top_neg <- gsea_df[order(gsea_df$NES), ][1:1, "ID"]

# Combine for plotting
combine <- c(top_pos, top_neg)

# Top 5 pathways
gseaplot2(gsea_results, geneSetID = 1:5, pvalue_table = TRUE)  ## use ful 

# Plot all 10 pathways in one overlay plot
top10<- gseaplot2(gsea_results,
                  geneSetID = combine,
                  pvalue_table = FALSE) 
top10 

ggsave(filename = "Classic_GSEA_male_vs_female.tiff",    
       plot = top10,                        
       width = 8, 
       height = 8,           
       dpi = 600,                          
       units = "in") 


# Gene Ontology 
library(clusterProfiler)
library(org.Hs.eg.db)

# ----- Gene Ranks -----
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Map SYMBOL -> ENTREZ
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# ----- Run GO GSEA -----
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

head(go_gsea@result)
gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_BP_male_vs_female.csv")

# ----- Add sign column (male = NES > 0, female = NES < 0) -----
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "male", "female")
go_gsea@result$sign <- go_df$sign

# ----- Dotplot split by NES sign -----
GO.dotplot <- dotplot(go_gsea, showCategory = 8, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Biological Process Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_male_vs_female_all.tiff",    
       plot = GO.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")

#Male old vs. Male young
#. Differential Expression using limma
# Design matrix
design <- model.matrix(~ 0 + factor(metadata$Old_young))
colnames(design) <- levels(factor(metadata$Old_young))

# Contrast male vs female
contrast.matrix <- makeContrasts(male_old_vs_female_old = m_old - f_old, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)

# Filter by adjusted p-value < 0.05
deg_filtered <- subset(results, adj.P.Val < 0.05)

# Check top rows
head(deg_filtered)

# Save to CSV
write.csv(deg_filtered, "DEG_male_vs_female_FDR0.05.csv", row.names = TRUE)

## Volcano plot 
library(ggplot2)
library(ggrepel)

# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
Male_old_vs_female_old<- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Downregulated" = "#00AFBB",
                                "Not significant" = "grey",
                                "Upregulated" = "#BB0C00")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Male_old vs Female_old") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                  size = 5, box.padding = 0.3, point.padding = 0.2,
                  max.overlaps = 20)

Male_old_vs_female_old

ggsave(filename = "Male_old_vs_female_old.tiff",    
       plot = Male_old_vs_female_old,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in") 

#compare old to young 

#. Differential Expression using limma
# Design matrix
design <- model.matrix(~ 0 + factor(metadata$Age_statutus))
colnames(design) <- levels(factor(metadata$Age_statutus))

# Contrast old vs. young
contrast.matrix <- makeContrasts(old_vs_young = old - young, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)

# Filter by adjusted p-value < 0.05
deg_filtered <- subset(results, adj.P.Val < 0.1)

# Check top rows
head(deg_filtered)

# Save to CSV
write.csv(deg_filtered, "DEG_male_vs_female_FDR0.05.csv", row.names = TRUE)

## Volcano plot 
library(ggplot2)
library(ggrepel)

# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
old_vs_young<- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Downregulated" = "#00AFBB",
                                "Not significant" = "grey",
                                "Upregulated" = "#BB0C00")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Old vs Young (all samples)") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                  size = 5, box.padding = 0.3, point.padding = 0.2,
                  max.overlaps = 20)
old_vs_young

ggsave(filename = "old_vs_young.tiff",    
       plot = old_vs_young,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in") 


# Design matrix
design <- model.matrix(~ 0 + factor(metadata$Old_young))
colnames(design) <- levels(factor(metadata$Old_young))

# Contrast male vs female
contrast.matrix <- makeContrasts(male_old_vs_male_old = m_old - m_young, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)

# Filter by adjusted p-value < 0.05
deg_filtered <- subset(results, adj.P.Val < 0.05)

# Check top rows
head(deg_filtered)

# Save to CSV
write.csv(deg_filtered, "DEG_male_vs_female_FDR0.05.csv", row.names = TRUE)

## Volcano plot 
library(ggplot2)
library(ggrepel)

# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
Male_old_vs_male_young<- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Downregulated" = "#00AFBB",
                                "Not significant" = "grey",
                                "Upregulated" = "#BB0C00")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Male_old vs male_young") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                  size = 5, box.padding = 0.3, point.padding = 0.2,
                  max.overlaps = 20)
Male_old_vs_male_young

ggsave(filename = "Male_old_vs_male_young.tiff",    
       plot = Male_old_vs_male_young,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in") 

## female old vs. female young
# Design matrix
design <- model.matrix(~ 0 + factor(metadata$Old_young))
colnames(design) <- levels(factor(metadata$Old_young))

# Contrast male vs female
contrast.matrix <- makeContrasts(female_old_vs_female_old = f_old - f_young, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)

# Filter by adjusted p-value < 0.05
deg_filtered <- subset(results, adj.P.Val < 0.05)

# Check top rows
head(deg_filtered)

# Save to CSV
write.csv(deg_filtered, "DEG_male_vs_female_FDR0.05.csv", row.names = TRUE)

## Volcano plot 
library(ggplot2)
library(ggrepel)

# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
female_old_vs_female_young<- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Downregulated" = "#00AFBB",
                                "Not significant" = "grey",
                                "Upregulated" = "#BB0C00")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("female_old vs female_young") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                  size = 5, box.padding = 0.3, point.padding = 0.2,
                  max.overlaps = 20)
female_old_vs_female_young


ggsave(filename = "female_old_vs_female_young.tiff",    
       plot = female_old_vs_female_young,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in") 


#male young vs. female young !! 
# Design matrix
design <- model.matrix(~ 0 + factor(metadata$Old_young))
colnames(design) <- levels(factor(metadata$Old_young))

# Contrast male vs female
contrast.matrix <- makeContrasts(male_young_vs_female_young = m_young - f_young, levels = design)

#fit the model
fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get all DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)

# Filter by adjusted p-value < 0.05
deg_filtered <- subset(results, adj.P.Val < 0.05)

# Check top rows
head(deg_filtered)

# Save to CSV
write.csv(deg_filtered, "DEG_male_vs_female_FDR0.05.csv", row.names = TRUE)

## Volcano plot 
library(ggplot2)
library(ggrepel)

# Use the full results table
all_genes <- results

# Add -log10 P-value for y-axis
all_genes$logP <- -log10(all_genes$P.Value)

# Add significance categories with fold-change cutoff
all_genes$Significance <- "Not significant"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC > 0.5] <- "Upregulated"
all_genes$Significance[all_genes$adj.P.Val < 0.05 & all_genes$logFC < -0.5] <- "Downregulated"

# Select top genes to label (e.g., top 10 up and top 10 down by adj.P.Val)
top_genes <- all_genes[all_genes$Significance != "Not significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 20)  # top 20 genes

# Volcano plot with gene labels
male_young_vs_female_young<- ggplot(all_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("Downregulated" = "#00AFBB",
                                "Not significant" = "grey",
                                "Upregulated" = "#BB0C00")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=13), 
        axis.title = element_text(face = "bold", size=13), 
        axis.text = element_text(face = "bold", size=13)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("male_young vs female_young") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                  size = 5, box.padding = 0.3, point.padding = 0.2,
                  max.overlaps = 20)
male_young_vs_female_young


ggsave(filename = "male_young_vs_female_young.tiff",    
       plot = male_young_vs_female_young,                        
       width = 8, height = 6,           
       dpi = 600,                          
       units = "in") 

#Enrichment analysis 

#Classic GSEA male vs.female 
# Create ranked list for GSEA (named vector: gene → score)
BiocManager::install("clusterProfiler")
gene_ranks <- results$t
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

library(clusterProfiler)
# download from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/ 
gmt_file <- "h.all.v7.5.1.symbols.gmt"
gene_sets <- read.gmt(gmt_file)

# Run GSEA
gsea_results <- GSEA(
  gene_ranks,
  TERM2GENE = gene_sets,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View results
head(gsea_results@result)
gsea_table <- as.data.frame(gsea_results)
write.csv(gsea_table, file = "hallmark_old_vs_young.csv")

# Add sign column==== very useful code 
hallmark_df <- as.data.frame(gsea_results)
hallmark_df$sign <- ifelse(hallmark_df$NES > 0, "old", "young")
gsea_results@result$sign <- hallmark_df$sign

# Dotplot split by NES sign
hallmark.dotplot<- dotplot(gsea_results, showCategory = 15, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hallmark Pathways")
hallmark.dotplot

#hallmark
ggsave(filename = "hallmark_old_vs_young_all.tiff",    
       plot = hallmark.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")

# ridge plot 
hallmark_male_vs_female<- ridgeplot(gsea_results) +
  theme_minimal() +
  theme(plot.title= element_text(hjust=0.5, face= "bold", size= 14),
        axis.text.y = element_text(size = 15, face= "bold"),  
        axis.text.x = element_text(size = 15, face= "bold"), 
        axis.title.x = element_text()) + 
  labs(x= "NES", title= "Male vs. Female")
hallmark_male_vs_female

#save plot
ggsave(filename = "hallmark_male_vs_female_all.tiff",    
       plot = hallmark_male_vs_female,                        
       width = 10, 
       height = 13,           
       dpi = 600,                          
       units = "in")

#visualization 
library(enrichplot)
library(ggplot2)

# Get top 5 enriched pathways (based on NES or p.adjust)
top_terms <- gsea_results@result$ID[1:5]

# Plot and save each one
for (i in seq_along(top_terms)) {
  p <- gseaplot2(gsea_results,
                 geneSetID = top_terms[i],
                 title = gsea_results@result$Description[i])
  
  # save as PDF or PNG
  ggsave(paste0("GSEA_plot_", i, ".pdf"), plot = p, width = 7, height = 5)
}

#or all result 
library(enrichplot)
# Convert GSEA results to data frame
gsea_df <- as.data.frame(gsea_results)
write.csv(gsea_df, file= "gsea_male_vs_female_allsample.csv")

# Top 5 positive NES
top_pos <- gsea_df[order(-gsea_df$NES), ][1:5, "ID"]

# Top 5 negative NES
top_neg <- gsea_df[order(gsea_df$NES), ][1:1, "ID"]

# Combine for plotting
combine <- c(top_pos, top_neg)

# Top 5 pathways
gseaplot2(gsea_results, geneSetID = 1:5, pvalue_table = TRUE)  ## use ful 

# Plot all 10 pathways in one overlay plot
top10<- gseaplot2(gsea_results,
                  geneSetID = combine,
                  pvalue_table = FALSE) 
top10 

ggsave(filename = "Classic_GSEA_male_vs_female.tiff",    
       plot = top10,                        
       width = 8, 
       height = 8,           
       dpi = 600,                          
       units = "in") 

#KEGG Pathway
install.packages("msigdbr") # if not installed
library(msigdbr)
library(GSEABase)

# Use t-statistics or logFC as ranking metric
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)


library(org.Hs.eg.db)
library(clusterProfiler)

# Map gene symbols to Entrez IDs
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# Keep only mapped genes
gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID

# Sort again just in case
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

## Run KEGG GSEA 
kegg_gsea <- gseKEGG(
  geneList     = gene_ranks,
  organism     = "hsa",          # human
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# View results
head(kegg_gsea@result)
gsea_table <- as.data.frame(kegg_gsea)
write.csv(gsea_table, file = "KEGG_old_vs_young.csv")

# Add sign column==== very useful code 
kegg_df <- as.data.frame(kegg_gsea)
kegg_df$sign <- ifelse(kegg_df$NES > 0, "old", "young")
kegg_gsea@result$sign <- kegg_df$sign

# Dotplot split by NES sign
KEGG.dotplot<- dotplot(kegg_gsea, showCategory = 12, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("KEGG Pathways")
KEGG.dotplot

#save plot
ggsave(filename = "KEGG_old_vs_young_all.tiff",    
       plot = KEGG.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


# Bar plot left and right =============
# Convert GSEA object to a data frame
gsea_df <- as.data.frame(kegg_gsea)

# Convert GSEA object to data frame
gsea_df <- as.data.frame(kegg_gsea)

# Top 10 positive NES
top_pos <- gsea_df[gsea_df$NES > 0, ]
top_pos <- top_pos[order(-top_pos$NES), ]
top_pos <- head(top_pos, 10)

# Top 10 negative NES
top_neg <- gsea_df[gsea_df$NES < 0, ]
top_neg <- top_neg[order(top_neg$NES), ]
top_neg <- head(top_neg, 10)

# Combine
top20 <- rbind(top_pos, top_neg)
top20$Description <- factor(top20$Description, levels = top20$Description[order(top20$NES)])

# Order pathways by NES for plotting  == new visualization 
ggplot(top20, aes(x = NES, y = Description, fill = NES > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#BB0C00", "FALSE" = "#00AFBB"),
                    labels = c("Negative NES", "Positive NES")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face="bold", size=12),
    axis.text.x = element_text(face="bold", size=12),
    plot.title = element_text(face="bold", size=14, hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_blank()
  ) + 
  theme(legend.position = "none") +
  labs(title = "Top 10 Positive & Top 10 Negative KEGG Pathways", x="NES", y="")

#save plot
ggsave(filename = "hallmark_male_vs_female_all.tiff",    
       plot = hallmark_male_vs_female,                        
       width = 10, 
       height = 13,           
       dpi = 600,                          
       units = "in")

#visualization 
library(enrichplot)
library(ggplot2)

# Get top 5 enriched pathways (based on NES or p.adjust)
top_terms <- gsea_results@result$ID[1:5]

# Plot and save each one
for (i in seq_along(top_terms)) {
  p <- gseaplot2(gsea_results,
                 geneSetID = top_terms[i],
                 title = gsea_results@result$Description[i])
  
  # save as PDF or PNG
  ggsave(paste0("GSEA_plot_", i, ".pdf"), plot = p, width = 7, height = 5)
}

#or all result 
library(enrichplot)
# Convert GSEA results to data frame
gsea_df <- as.data.frame(gsea_results)
write.csv(gsea_df, file= "gsea_male_vs_female_allsample.csv")

# Top 5 positive NES
top_pos <- gsea_df[order(-gsea_df$NES), ][1:5, "ID"]

# Top 5 negative NES
top_neg <- gsea_df[order(gsea_df$NES), ][1:1, "ID"]

# Combine for plotting
combine <- c(top_pos, top_neg)

# Top 5 pathways
gseaplot2(gsea_results, geneSetID = 1:5, pvalue_table = TRUE)  ## use ful 

# Plot all 10 pathways in one overlay plot
top10<- gseaplot2(gsea_results,
                  geneSetID = combine,
                  pvalue_table = FALSE) 
top10 

ggsave(filename = "Classic_GSEA_male_vs_female.tiff",    
       plot = top10,                        
       width = 8, 
       height = 8,           
       dpi = 600,                          
       units = "in") 


# Gene Ontology 
library(clusterProfiler)
library(org.Hs.eg.db)

# ----- Gene Ranks -----
gene_ranks <- results$t             # or results$logFC
names(gene_ranks) <- rownames(results)

# Map SYMBOL -> ENTREZ
gene_mapping <- bitr(names(gene_ranks), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

gene_ranks <- gene_ranks[gene_mapping$SYMBOL]
names(gene_ranks) <- gene_mapping$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

#  Run GO GSEA -----
go_gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",        # can be "BP", "MF", or "CC"
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

head(go_gsea@result)
gsea_GO_table <- as.data.frame(go_gsea)
write.csv(gsea_GO_table, file = "GO_BP_old_vs_young.csv")

# ----- Add sign column (Old = NES > 0, Young = NES < 0) -----
go_df <- as.data.frame(go_gsea)
go_df$sign <- ifelse(go_df$NES > 0, "Old", "Young")
go_gsea@result$sign <- go_df$sign

# ----- Dotplot split by NES sign -----
GO.dotplot <- dotplot(go_gsea, showCategory = 7, split = "sign") +
  facet_grid(. ~ sign) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Biological Process Pathways")

GO.dotplot

#save plot
ggsave(filename = "GO_old_vs_young_all.tiff",    
       plot = GO.dotplot,                        
       width = 8, 
       height = 9,           
       dpi = 600,                          
       units = "in")


### MCP counter 
BiocManager::install("MCPcounter")
require(MCPcounter)
install.packages(c("e1071", "preprocessCore", "tidyverse"))
library(e1071)
library(preprocessCore)
library(tidyverse)

# Save CIBERSORT-ready file
write.table(exprs_gene, "exprs_CIBERSORT.txt", sep = "\t", quote = FALSE, col.names = NA)

#Run Cibersort 
# Load CIBERSORT script
source("CIBERSORT.R")

# Run CIBERSORT
cibersort_results <- CIBERSORT("LM22.txt", "exprs_CIBERSORT.txt", perm = 100, QN = TRUE)

# View results
head(cibersort_results)

# Assume exprs_gene is a matrix/dataframe with genes in rownames
exprs_to_upload <- data.frame(GeneSymbol = rownames(exprs_gene), exprs_gene)

# Write out with proper header
write.table(exprs_to_upload, 
            "exprs_CIBERSORT.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)


#' read cibersort data ===============================================================================================
cibersort<- read.csv("Cibersort.csv", 
                     stringsAsFactors = TRUE)


#using pheatmap 
require(tidyheatmaps)

#prepare 
df_long<- cibersort %>% 
  pivot_longer(
    cols = naive_B_cell:neutrophils, 
    names_to = "cell_type", 
    values_to = "value"
  )

## heatmap 
tidyheatmap(df_long, 
            rows = cell_type, 
            columns = Tumor,
            values = value,
            scale = "row",
            annotation_col = c(Gender),
            gaps_col = Gender)

df_long$new <- df_long$value * 100

#stack bar graph for all sample 
library(ggplot2)

# Assign specific colors to some cell types
cibersort_colors <- c(
  "#FF0000","#FF6666","#FF9999",
  "#0000FF","#6666FF","#9999FF","#000000",
  "#00CCFF","#00FFFF","#339999",
  "#00FF00","#66FF66",
  "#CC9900","#996600","#FFCC00","#FF9900",
  "#6600CC","#9933FF",
  "#FFCCFF","#FF99FF",
  "#FF9966","#FF3300"
)

df_long$cell_type<- factor(df_long$cell_type, 
                           levels = c("naive_B_cell",
                                      "memory_B_cell",
                                      "plasma_cell",
                                      "CD8",
                                      "naive_CD4",
                                      "resting_CD4_memory",
                                      "activated_memory_CD4",
                                      "Tfh",
                                      "Treg",
                                      "gdT",
                                      "resting_NK",
                                      "activated_NK",
                                      "monocytes",
                                      "M0_Macrophages",
                                      "M1_Macrophages",
                                      "M2_Macrophages",
                                      "resting_DC",
                                      "activated_DC",
                                      "resting_Mast_cell",
                                      "activated_Mast_cell",
                                      "eosinophils",
                                      "neutrophils"))

df_long$Tumor<- factor(df_long$Tumor, 
                       levels = c("ATC_1", 
                                  "ATC_2", 
                                  "ATC_3", 
                                  "ATC_4", 
                                  "ATC_5",
                                  "ATC_6", 
                                  "ATC_7", 
                                  "ATC_8", 
                                  "ATC_9", 
                                  "ATC_10",
                                  "ATC_11", 
                                  "ATC_12", 
                                  "ATC_13", 
                                  "ATC_14", 
                                  "ATC_15",
                                  "ATC_16", 
                                  "ATC_17", 
                                  "ATC_18", 
                                  "ATC_19", 
                                  "ATC_20"))
## Plot
fraction<- ggplot(df_long, aes(x = Old_young_n,
                               y = new, 
                               fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cibersort_colors) +
  theme_classic() +
  labs(
    title = "Fraction",
    x = "Tumor",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold", 
                                    size= 12), 
        axis.text = element_text(size= 12, 
                                 face = "bold"))
fraction

ggsave("fraction.ATC.tiff", 
       plot = fraction,
       height=7, 
       width = 13, 
       units = "in", 
       dpi = 300)


#male vs. female 
library(dplyr)
library(ggplot2)

# Compute average fraction per Gender and cell_type
df_avg <- df_long %>%
  group_by(Gender, cell_type) %>%
  summarise(avg_value = mean(new), .groups = "drop")  # 'new' is value * 100

# Stacked bar plot with averages
fraction_m_vs_f_60<- ggplot(df_avg, aes(x = Gender, 
                                         y = avg_value, 
                                         fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cibersort_colors) +
  theme_classic() +
  labs(
    title = "Average Fraction by Gender",
    x = "Gender",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold")
  )
fraction_m_vs_f_60

ggsave("fraction_m_vs_f_60.tiff", 
       plot = fraction_m_vs_f_60,
       height=7, 
       width = 7, 
       units = "in", 
       dpi = 300)

# boxplot for invidividual cells 
boxplot_m_vs_f<- ggplot(df_long, aes(x= cell_type, 
                    y= new, 
                    fill = Gender)) + 
  geom_boxplot() + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, 
                               hjust = 1, 
                               size= 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold"), 
    legend.title = element_blank()) + 
  labs(y= "% of CD45", title= "male vs. female")
boxplot_m_vs_f

ggsave("fraction_m_vs_f_individual.tiff", 
       plot = boxplot_m_vs_f,
       height=4, 
       width = 7, 
       units = "in", 
       dpi = 300)

#older vs_young 

# Compute average fraction per Gender and cell_type
df_avg <- df_long %>%
  group_by(Age_statutus, cell_type) %>%
  summarise(avg_value = mean(new), .groups = "drop")  # 'new' is value * 100

# Stacked bar plot with averages
fraction_old_vs_young_60 <- ggplot(df_avg, aes(x = Age_statutus, y = avg_value, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cibersort_colors) +
  theme_classic() +
  labs(
    title = "Average Fraction by Age",
    x = "Age",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold")
  )

fraction_old_vs_young_60

ggsave("fraction_old_vs_young_60.tiff", 
       plot = fraction_old_vs_young_60,
       height=7, 
       width = 7, 
       units = "in", 
       dpi = 300)

boxplot_old_vs_young_60<- ggplot(df_long, aes(x= cell_type, 
                                     y= new, 
                                     fill = Age_statutus)) + 
  geom_boxplot() + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, 
                               hjust = 1, 
                               size= 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold"), 
    legend.title = element_blank()) + 
  labs(y= "% of CD45", title= "Old vs. young")
boxplot_old_vs_young_60

ggsave("fraction_old_vs_young_individual.tiff", 
       plot = boxplot_old_vs_young_60,
       height=4, 
       width = 7, 
       units = "in", 
       dpi = 300)


## Male_old vs. Male_young
Male_old_vs_male_young<- df_long %>% filter(Old_young %in% c("m_old", 
                                                             "m_young"))


df_avg <- Male_old_vs_male_young %>%
  group_by(Old_young, cell_type) %>%
  summarise(avg_value = mean(new), .groups = "drop")  # 'new' is value * 100

# Stacked bar plot with averages
fraction_male_old_vs_male_young <- ggplot(df_avg, aes(x = Old_young, 
                                                      y = avg_value, 
                                                      fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cibersort_colors) +
  theme_classic() +
  labs(
    title = "Old vs. young (male)",
    x = "Age",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5 ),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold")
  )

fraction_male_old_vs_male_young


ggsave("fraction_male_old_vs_male_young.tiff", 
       plot = fraction_male_old_vs_male_young,
       height=7, 
       width = 7, 
       units = "in", 
       dpi = 300)


old_vs_young_male<- ggplot(Male_old_vs_male_young, aes(x= cell_type, 
                    y= new, 
                    fill = Old_young)) + 
  geom_boxplot() + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, 
                               hjust = 1, 
                               size= 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold"), 
    legend.title = element_blank()) + 
  labs(y= "% of CD45", title= "Old vs. young (male)")
old_vs_young_male

ggsave("old_vs_young_male_individual.tiff", 
       plot = old_vs_young_male,
       height=4, 
       width = 7, 
       units = "in", 
       dpi = 300)


## female_old vs. female_young
female_old_vs_female_young<- df_long %>% filter(Old_young %in% c("f_old", 
                                                             "f_young"))


df_avg <- female_old_vs_female_young %>%
  group_by(Old_young, cell_type) %>%
  summarise(avg_value = mean(new), .groups = "drop")  # 'new' is value * 100

# Stacked bar plot with averages
fraction_female_old_vs_female_young <- ggplot(df_avg, aes(x = Old_young, 
                                                      y = avg_value, 
                                                      fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cibersort_colors) +
  theme_classic() +
  labs(
    title = "Old vs. young (female)",
    x = "Age",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5 ),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold")
  )
fraction_female_old_vs_female_young 


ggsave("fraction_female_old_vs_female_young.tiff", 
       plot = fraction_female_old_vs_female_young,
       height=7, 
       width = 7, 
       units = "in", 
       dpi = 300)


old_vs_young_female<- ggplot(female_old_vs_female_young, aes(x= cell_type, 
                                                       y= new, 
                                                       fill = Old_young)) + 
  geom_boxplot() + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, 
                               hjust = 1, 
                               size= 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold"), 
    legend.title = element_blank()) + 
  labs(y= "% of CD45", title= "Old vs. young (female)")
old_vs_young_female

ggsave("old_vs_young_female_individual.tiff", 
       plot = old_vs_young_female,
       height=4, 
       width = 7, 
       units = "in", 
       dpi = 300)


## Male old vs_ female old 
male_old_vs_female_old<- df_long %>% filter(Old_young %in% c("m_old", 
                                                                 "f_old"))


df_avg <- male_old_vs_female_old %>%
  group_by(Old_young, cell_type) %>%
  summarise(avg_value = mean(new), .groups = "drop")  # 'new' is value * 100

# Stacked bar plot with averages
fraction_male_old_vs_female_old <- ggplot(df_avg, aes(x = Old_young, 
                                                          y = avg_value, 
                                                          fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cibersort_colors) +
  theme_classic() +
  labs(
    title = "male_Old vs. female_old",
    x = "Age",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5 ),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold")
  )
fraction_male_old_vs_female_old


ggsave("fraction_male_old_vs_female_old.tiff", 
       plot = fraction_male_old_vs_female_old,
       height=7, 
       width = 7, 
       units = "in", 
       dpi = 300)


male_old_vs_female_old<- ggplot(male_old_vs_female_old, aes(x= cell_type, 
                                                             y= new, 
                                                             fill = Old_young)) + 
  geom_boxplot() + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, 
                               hjust = 1, 
                               size= 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold"), 
    legend.title = element_blank()) + 
  labs(y= "% of CD45", title= "Male_old vs. female_old")
male_old_vs_female_old

ggsave("male_old_vs_female_old_individual.tiff", 
       plot = male_old_vs_female_old,
       height=4, 
       width = 7, 
       units = "in", 
       dpi = 300)

# male_young_vs.female_young
male_young_vs_female_young<- df_long %>% filter(Old_young %in% c("m_young", 
                                                             "f_young"))


df_avg <- male_young_vs_female_young %>%
  group_by(Old_young, cell_type) %>%
  summarise(avg_value = mean(new), .groups = "drop")  # 'new' is value * 100

# Stacked bar plot with averages
fraction_male_young_vs_female_young <- ggplot(df_avg, aes(x = Old_young, 
                                                      y = avg_value, 
                                                      fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cibersort_colors) +
  theme_classic() +
  labs(
    title = "male_young vs. female_young",
    x = "Age",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5 ),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold")
  )
fraction_male_young_vs_female_young


ggsave("fraction_male_young_vs_female_young.tiff", 
       plot = fraction_male_young_vs_female_young,
       height=7, 
       width = 7, 
       units = "in", 
       dpi = 300)


male_young_vs_female_young<- ggplot(male_young_vs_female_young, 
                                    aes(x= cell_type, 
                                    y= new, 
                                    fill = Old_young)) + 
  geom_boxplot() + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, 
                               hjust = 1, 
                               size= 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size=12),
    axis.text = element_text(size=12, face="bold"), 
    legend.title = element_blank()) + 
  labs(y= "% of CD45", title= "male_young vs. female_young")
male_young_vs_female_young

ggsave("male_young_vs_female_young_individual.tiff", 
       plot = male_young_vs_female_young,
       height=4, 
       width = 7, 
       units = "in", 
       dpi = 300)

#perform classic GSEA 
# Create ranked list for GSEA (named vector: gene score)
gene_ranks <- deg_results$t
names(gene_ranks) <- rownames(deg_results)

# Sort decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
head(gene_ranks)

#' enf of Microarray data analysis for anaplatic thyroid cancer