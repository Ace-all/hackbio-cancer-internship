\# Load required libraries

library(TCGAbiolinks)

library(SummarizedExperiment)

library(DESeq2)

library(clusterProfiler)

library(org.Hs.eg.db)

library(randomForest)

library(ggplot2)

library(dplyr)

library(pheatmap)

\#Download TCGA-BRCA dataset

query <- GDCquery(

  project = "TCGA-BRCA",

  data.category = "Transcriptome Profiling",

  data.type = "Gene Expression Quantification",

  workflow\.type = "STAR - Counts",

                  access = "open")

GDCdownload(query)

data <- GDCprepare(query)

\#Clean and preprocess the data

\# Extract counts,metadata and sample information

counts <- assay(data)

metadata <- as.data.frame(colData(data))

sample\_info <- as.data.frame(colData(data))

\# Handle missing values and normalize data

dds <- DESeqDataSetFromMatrix(countData = counts,

                              colData = sample\_info,

                              design = \~ 1)

dds <- estimateSizeFactors(dds)

normalized\_counts <- counts(dds, normalized = TRUE)

\# Select samples based on available sample types

sample\_types <- unique(sample\_info$sample\_type)

print("Available sample types:")

print(table(sample\_info$sample\_type))

\# Function to select samples

select\_samples <- function(sample\_info, sample\_type, n = 20) {

  samples <- which(sample\_info$sample\_type == sample\_type)

  if (length(samples) > n) {

    samples <- sample(samples, n)

  }

  return(samples)

}

\# Select samples (adjust these based on available sample types)

primary\_samples <- select\_samples(sample\_info, "Primary Tumor")

solid\_tissue\_normal\_samples <- select\_samples(sample\_info, "Solid Tissue Normal")

selected\_samples <- c(primary\_samples, solid\_tissue\_normal\_samples)

selected\_counts <- normalized\_counts\[, selected\_samples]

selected\_info <- sample\_info\[selected\_samples, ]

\# Differential expression analysis

dds\_selected <- DESeqDataSetFromMatrix(countData = round(selected\_counts),

                                       colData = selected\_info,

                                       design = \~ sample\_type)

dds\_selected <- DESeq(dds\_selected)

res <- results(dds\_selected)

\# Sort results by adjusted p-value

res\_ordered <- res\[order(res$padj), ]

sig\_genes <- subset(res\_ordered, padj < 0.05 & abs(log2FoldChange) > 1)

\# Functional enrichment analysis

entrez\_ids <- mapIds(org.Hs.eg.db, keys = rownames(sig\_genes), keytype = "ENSEMBL", column = "ENTREZID")

  # GO Enrichment Analysis

  go\_enrichment <- enrichGO(gene = entrez\_ids,

                            OrgDb = org.Hs.eg.db,

                            ont = "BP",

                            pAdjustMethod = "BH",

                            pvalueCutoff = 0.05,

                            qvalueCutoff = 0.05)

  

  # Print GO enrichment results

  cat("\nTop 10 enriched GO terms:\n")

  print(head(summary(go\_enrichment), 10))

  

  # Save GO enrichment results

  write.csv(summary(go\_enrichment), "go\_enrichment\_results.csv")

  

  # KEGG Pathway Enrichment Analysis

  kegg\_enrichment <- enrichKEGG(gene = entrez\_ids,

                                organism = "hsa",

                                pvalueCutoff = 0.05)

  

  # Print KEGG enrichment results

  cat("\nTop 10 enriched KEGG pathways:\n")

  print(head(summary(kegg\_enrichment), 10))

  

  # Save KEGG enrichment results

  write.csv(summary(kegg\_enrichment), "kegg\_enrichment\_results.csv")

} else {

  cat("No valid ENTREZ IDs found. Skipping enrichment analysis.\n")

}

entrez\_ids <- mapIds(org.Hs.eg.db, keys = rownames(sig\_genes), keytype = "ENSEMBL", column = "ENTREZID")

go\_enrichment <- enrichGO(gene = entrez\_ids,

                          OrgDb = org.Hs.eg.db,

                          ont = "BP",

                          pAdjustMethod = "BH",

                          pvalueCutoff = 0.05,

                          qvalueCutoff = 0.05)

\#Machine Learning

\# Prepare data for ML

all.trans <- data.frame(t(selected\_counts))

\# Select the top 1000 most variable genes based on standard deviation

SDs <- apply(all.trans, 2, sd)

topPreds <- order(SDs, decreasing = TRUE)\[1:1000]

all.trans <- all.trans\[, topPreds]

\# PreProcessing steps

\# Remove near-zero variance predictors

all.zero <- preProcess(all.trans, method = 'nzv', uniqueCut = 15)

all.trans <- predict(all.zero, all.trans)

\# Center the data

all.center <- preProcess(all.trans, method = 'center')

all.trans <- predict(all.center, all.trans)

\# Remove highly correlated features

all.corr <- preProcess(all.trans, method = 'corr', cutoff = 0.5)

all.trans <- predict(all.corr, all.trans)

ml\_data <- all.trans

ml\_labels <- factor(selected\_info$sample\_type)

\# Feature selection (using variance)

var\_genes <- apply(ml\_data, 2, var)

top\_genes <- names(sort(var\_genes, decreasing = TRUE))\[1:1000]

ml\_data\_filtered <- ml\_data\[, top\_genes]

\# Split data into training and testing sets

set.seed(42)

train\_indices <- sample(1:nrow(ml\_data\_filtered), 0.7 \* nrow(ml\_data\_filtered))

train\_data <- ml\_data\_filtered\[train\_indices, ]

train\_labels <- ml\_labels\[train\_indices]

test\_data <- ml\_data\_filtered\[-train\_indices, ]

test\_labels <- ml\_labels\[-train\_indices]

\# Random Forest classification

rf\_model <- randomForest(x = train\_data, y = train\_labels, ntree = 500)

rf\_predictions <- predict(rf\_model, test\_data)

rf\_accuracy <- sum(rf\_predictions == test\_labels) / length(test\_labels)

conf\_matrix <- confusionMatrix(rf\_predictions, test\_labels)

conf\_matrix\_plot <- as.data.frame(conf\_matrix$table)

conf\_matrix\_plot$Prediction <- factor(conf\_matrix\_plot$Prediction, levels = rev(levels(conf\_matrix\_plot$Prediction)))

\# Visualizations

\# Volcano plot

volcano\_plot <- ggplot(data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +

  geom\_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.6) +

  scale\_color\_manual(values = c("grey", "red")) +

  labs(title = "Volcano Plot of Differential Expression",

       x = "log2 Fold Change",

       y = "-log10 Adjusted p-value") +

  theme\_minimal()

\# confusion matrix

ggplot(conf\_matrix\_plot, aes(Prediction, Reference, fill = Freq)) +

  geom\_tile() +

  geom\_text(aes(label = sprintf("%d", Freq)), vjust = 1) +

  scale\_fill\_gradient(low = "white", high = "steelblue") +

  theme\_minimal() +

  labs(title = "Confusion Matrix",

       x = "Predicted",

       y = "Actual")

\# ROC Curve

roc\_obj <- roc(test\_labels, as.numeric(rf\_predictions))

roc\_plot <- ggroc(roc\_obj) +

  geom\_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray") +

  theme\_minimal() +

  labs(title = paste("ROC Curve (AUC =", round(auc(roc\_obj), 2), ")"),

       x = "False Positive Rate",

       y = "True Positive Rate")

ggsave("roc\_curve\_plot.png", roc\_plot, width = 8, height = 6)

\# Feature Importance

importance\_scores <- importance(rf\_model)

top\_features <- head(importance\_scores\[order(importance\_scores\[, 1], decreasing = TRUE), , drop = FALSE], 20)

feature\_importance\_plot <- ggplot(data.frame(feature = rownames(top\_features), importance = top\_features\[, 1]),

                                  aes(x = reorder(feature, importance), y = importance)) +

  geom\_bar(stat = "identity", fill = "steelblue") +

  coord\_flip() +

  theme\_minimal() +

  labs(title = "Top 20 Feature Importance",

       x = "Features",

       y = "Importance")

ggsave("feature\_importance\_plot.png", feature\_importance\_plot, width = 10, height = 8)

\# Heatmap of top differentially expressed genes

top\_genes <- head(rownames(res\_ordered), 50)

heatmap\_data <- selected\_counts\[top\_genes, ]

heatmap\_plot <- pheatmap(heatmap\_data,

                         scale = "row",

                         annotation\_col = data.frame(sample\_type = selected\_info$sample\_type),

                         show\_rownames = FALSE,

                         main = "Heatmap of Top 50 Differentially Expressed Genes")

\# Save results and plots

write.csv(res\_ordered, "differential\_expression\_results.csv")

write.csv(summary(go\_enrichment), "go\_enrichment\_results.csv")

ggsave("volcano\_plot.png", volcano\_plot)



dev.off()

\# Print summary of results

cat("Number of significant differentially expressed genes:", nrow(sig\_genes), "\n")

cat("Top 5 enriched GO terms:\n")

print(head(summary(go\_enrichment), 5))

cat("Random Forest classification accuracy:", rf\_accuracy, "\n")
