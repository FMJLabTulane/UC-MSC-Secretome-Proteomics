# Load necessary libraries
library(tidyverse)
library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)


# 2. DATA LOADING ----
# Load cleaned data and metadata
proteomics_metadata <- read.csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/Solo analysis using Cade only Hu data/1_proteomics_metadata.csv")
metadata <- read.csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/samples_metadata.csv")


#3. Separate multiple protein and IBAQ values in a row

# Step 1: Identify metadata columns (those before the intensity columns)
meta_cols <- grep("PG\\.", names(proteomics_metadata)
                  
                  # Step 2: Separate the semicolon-delimited proteins in "PG.ProteinGroups"
                  
                  # Step 1: Add row ID
                  proteomics_metadata <- proteomics_metadata %>%
                    mutate(row_id = row_number())
                  
                  # Step 2: Separate annotation columns (expanding rows by semicolons)
                  annotation_cols <- c("PG.ProteinNames", "PG.ProteinGroups", "PG.Genes", "PG.Organisms", "PG.ProteinDescriptions","PG.UniProtIds")
                  
                  proteomics_long <- proteomics_metadata %>%
                    separate_rows(all_of(annotation_cols), sep = ";") %>%
                    group_by(row_id) %>%
                    mutate(protein_index = row_number()) %>%
                    ungroup()
                  
                  # Step 3: Expand IBAQ values across samples
                  ibaq_cols <- grep("IBAQ", names(proteomics_metadata), value = TRUE)
                  
                  ibaq_long <- proteomics_metadata %>%
                    select(row_id, all_of(ibaq_cols)) %>%
                    pivot_longer(-row_id, names_to = "Sample_ID", values_to = "IBAQ_raw") %>%
                    mutate(Sample_ID_clean = str_remove(Sample_ID, "\\.raw\\.PG\\.IBAQ$")) %>%
                    separate_rows(IBAQ_raw, sep = ";") %>%
                    group_by(row_id, Sample_ID_clean) %>%
                    mutate(protein_index = row_number()) %>%
                    ungroup() %>%
                    mutate(IBAQ = as.numeric(IBAQ_raw)) %>%
                    select(row_id, protein_index, Sample_ID_clean, IBAQ)
                  
                  # Step 4: Join annotations with IBAQ values
                  merged_long <- proteomics_long %>%
                    left_join(ibaq_long, by = c("row_id", "protein_index"))
                  
                  # Step 5: Convert to wide format: One row per protein, one column per sample
                  merged_wide <- merged_long %>%
                    select(PG.ProteinGroups, PG.Genes, PG.ProteinDescriptions, PG.UniProtIds, PG.ProteinNames,
                           Sample_ID_clean, IBAQ) %>%
                    pivot_wider(names_from = Sample_ID_clean, values_from = IBAQ)
                  
                  # View result
                  head(merged_wide)
                  
                  # Save the final cleaned dataframe
                  write.csv(merged_wide, "cleaned_proteomics.csv", row.names = FALSE)
                  
                  # 2. DATA LOADING ----
                  # Load cleaned data and metadata
                  cleaned_proteomics <- read.csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/cleaned_proteomics.csv")
                  metadata <- read.csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/fq_proteomics_analysis/rawdata/iBAQ/samples_metadata.csv")
                  
                  head(cleaned_proteomics, n=5)
                  
                  
                  
                  # Remove annotation columns to calculate missing per row
                  data_numeric <- cleaned_proteomics[, 6:ncol(cleaned_proteomics)]
                  
                  # Step 2: Filter proteins with too many NAs (e.g. >50%)
                  keep_rows <- rowMeans(is.na(data_numeric)) <= 0.5
                  filtered_proteomics <- cleaned_proteomics[keep_rows, ]
                  
                  # Step 3: Optional - Impute NAs with half-min of each column (simple imputation)
                  imputed_data <- filtered_proteomics
                  imputed_data[, 6:ncol(imputed_data)] <- lapply(imputed_data[, 6:ncol(imputed_data)], function(x) {
                    x[is.na(x)] <- min(x, na.rm = TRUE) / 2
                    return(x)
                  })
                  
                  # Step 4: Save cleaned version
                  write.csv(imputed_data, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/imputed_proteomics.csv", row.names = FALSE)
                  
                  imputed_data <- read.csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/Solo analysis using Cade only Hu data/imputed_proteomics.csv")
                  
                  
                  
                  # 1. LOAD CLEANED & IMPUTED DATA ----
                  log2_proteomics <- imputed_data  # <- make sure this is loaded in your environment
                  log2_proteomics[, 6:ncol(log2_proteomics)] <- log2(log2_proteomics[, 6:ncol(log2_proteomics)] + 0.1)
                  
                  # Step 4: Save Log2 version
                  write.csv(log2_proteomics, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/log2_proteomics.csv", row.names = FALSE)
                  
                  # Step 4: Save cleaned version
                  write.csv(imputed_data, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/imputed_proteomics.csv", row.names = FALSE)
                  
                  
                  
                  
                  
########### NH have edited the log2_proteomics data according to the FQ cleaned proteomics data: column number and the column orientation were not same and needs to be aligned for downstream analysis####

log2_proteomics <- read.csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/Solo analysis using Cade only Hu data/log2_proteomics.csv")
metadata <- read.csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/Solo analysis using Cade only Hu data/samples_metadata.csv")
                  

# Create expression matrix
expr_matrix <- log2_proteomics %>% 
  select(-(1:6)) %>% 
  as.matrix()
rownames(expr_matrix) <- log2_proteomics$PG.ProteinNames

# 3. QUALITY CONTROL ----
# PCA plot to check sample clustering
pca <- prcomp(t(expr_matrix), scale. = TRUE)
pca_data <- as.data.frame(pca$x) %>% 
  mutate(Sample = rownames(.)) %>% 
  left_join(metadata, by = c("Sample" = "SampleID"))

ggplot(pca_data, aes(PC1, PC2, color = Group, shape = Sex)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Protein Expression",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))

# PCA plot to check sample clustering for each sex separately
ggplot(pca_data, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  facet_wrap(~ Sex) +  # Creates separate panels
  theme_minimal() +
  labs(title = "PCA of Protein Expression by Sex",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  theme(strip.text = element_text(size = 12))  # Larger facet labels

# Removing outliers from PCA (if you have):Calculate Mahalanobis distance (multivariate outlier detection)
pca_scores <- pca_data[, c("PC1", "PC2")]
mahalanobis_dist <- mahalanobis(pca_scores, 
                                center = colMeans(pca_scores), 
                                cov = cov(pca_scores))

# Identify outliers (e.g., beyond 97.5% quantile of Chi-squared distribution)
threshold <- qchisq(0.975, df = 2)  # df = number of PCs considered
outliers <- which(mahalanobis_dist > threshold)

#Visualization of outliers
ggplot(pca_data, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  geom_point(data = pca_data[outliers, ], 
             aes(PC1, PC2), color = "black", size = 5, shape = 1) +
  labs(title = "PCA with Outliers Highlighted")

# Remove outliers from expression matrix
expr_matrix_clean <- expr_matrix[-outliers, ]

# Remove corresponding metadata
proteomics_data_clean <- log2_proteomics[-outliers, ]
pca_data_clean <- pca_data[-outliers, ]

# Cleaned PCA plot after outlier removal for each sex separately
ggplot(pca_data_clean, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  facet_wrap(~ Sex) +  # Creates separate panels
  theme_minimal() +
  labs(title = "PCA of PE by Sex After Outlier Removal",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  theme(strip.text = element_text(size = 12))  # Larger facet labels





# 4. DIFFERENTIAL EXPRESSION ----
# Design matrix
groups <- factor(paste(metadata$Sex, metadata$Group, sep = "_"))
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

# Contrasts
contrasts <- makeContrasts(
  Male_DHT = Male_DHT - Male_EtOH,
  Female_DHT = Female_DHT - Female_EtOH,
  Male_T = Male_T - Male_EtOH,
  Female_T = Female_T - Female_EtOH,
  Male_E2 = Male_E2 - Male_EtOH,
  Female_E2 = Female_E2 - Female_EtOH,
  Male_TA = Male_T_A - Male_EtOH,
  Female_TA = Female_T_A - Female_EtOH,
  Male_TS = Male_T_S - Male_EtOH,
  Female_TS = Female_T_S - Female_EtOH,
  levels = design
)

# Fit model
fit <- lmFit(expr_matrix_clean, design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# 5. VISUALIZATION ----
## Volcano plots for each comparison
generate_volcano <- function(comparison, fit_object) {
  results <- topTable(fit_object, coef = comparison, number = Inf)
  
  EnhancedVolcano(results,
                  lab = rownames(results),
                  x = 'logFC',
                  y = 'P.Value',
                  title = comparison,
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 2.0,
                  labSize = 4.0,
                  colAlpha = 0.6)
}

# Generate all volcano plots
volcano_plots <- map(colnames(contrasts), ~generate_volcano(., fit2))
names(volcano_plots) <- colnames(contrasts)

# Save volcano plots
walk2(volcano_plots, names(volcano_plots), 
      ~ggsave(paste0("volcano_", .y, ".png"), .x, width = 8, height = 8))

getwd()

## Heatmap of top differentially expressed proteins
# Get top proteins from any comparison
top_proteins <- topTable(fit2, number = 50, sort.by = "B")$ID

# Z-score normalize expression
heatmap_data <- expr_matrix_clean[top_proteins, ] %>% 
  t() %>% 
  scale() %>% 
  t()

# Annotation
ha <- HeatmapAnnotation(
  Group = metadata$Group,
  Sex = metadata$Sex,
  col = list(
    Group = c("EtOH" = "grey", "DHT" = "blue", "E2" = "pink", 
              "T" = "orange", "T_A" = "red", "T_S" = "purple"),
    Sex = c("Male" = "darkblue", "Female" = "deeppink")
  )
)

# Plot heatmap
Heatmap(heatmap_data,
        name = "Z-score",
        top_annotation = ha,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 8),
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

# 6. IMPROVED PROTEIN EXPRESSION PLOTS ----
plot_protein <- function(protein, data = proteomics_data_clean, metadata = metadata) {
  # Extract expression values
  expr <- data %>% 
    filter(PG.ProteinNames == protein) %>% 
    select(matches("^m_|^f_")) %>% 
    gather(key = "Sample", value = "Expression")
  
  # Add metadata
  expr <- expr %>% 
    left_join(metadata, by = c("Sample" = "SampleID")) %>% 
    mutate(Group = factor(Group, levels = c("EtOH", "E2", "DHT", "T", "T_A", "T_S")),
           Sex = factor(Sex, levels = c("Male", "Female")))
  
  # Plot
  ggplot(expr, aes(x = Group, y = Expression, fill = Sex)) +
    geom_boxplot(position = position_dodge(0.8), alpha = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), 
               size = 2, alpha = 0.7) +
    scale_fill_manual(values = c("Male" = "blue", "Female" = "pink")) +
    labs(title = protein,
         y = "log2(iBAQ)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Example usage:
# plot_protein("GFRA2")

# 7. SAVE ALL RESULTS ----
# Function to save complete DPE results
save_dpe_results <- function(comparison, fit_object) {
  results <- topTable(fit_object, coef = comparison, number = Inf)
  
  # Add protein annotations
  results <- results %>% 
    rownames_to_column("Protein") %>% 
    left_join(proteomics_data_clean %>% 
                select(PG.ProteinNames, PG.Genes, PG.ProteinDescriptions) %>% 
                distinct(),
              by = c("Protein" = "PG.ProteinNames"))
  
  write.csv(results, paste0("DPE_results_", comparison, ".csv"), row.names = FALSE)
}

# Save all comparisons
walk(colnames(contrasts), ~save_dpe_results(., fit2))

####### 8. CLUSTERED HEATMAP (SEX & GROUP) ----
# Prepare expression matrix (use cleaned data)
heatmap_data <- expr_matrix_clean

# Z-score normalization (row-wise)
heatmap_z <- t(scale(t(heatmap_data)))

# Annotation (ensure metadata matches expr_matrix columns)
ha <- HeatmapAnnotation(
  Group = metadata$Group,
  Sex = metadata$Sex,
  col = list(
    Group = c("EtOH" = "grey", "DHT" = "blue", "E2" = "pink", 
              "T" = "orange", "T_A" = "red", "T_S" = "purple"),
    Sex = c("Male" = "darkblue", "Female" = "deeppink")
  ),
  annotation_name_side = "left"  # Moves annotation titles to the left
)

# Plot heatmap with adjusted legend
Heatmap(heatmap_z,
        name = "Z-score",
        top_annotation = ha,
        column_split = metadata$Sex,
        show_column_names = FALSE,
        show_row_names = FALSE,
        cluster_columns = TRUE,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        heatmap_legend_param = list(
          title_position = "topcenter",    # Valid position for horizontal legend
          legend_direction = "vertical",
          legend_width = unit(5, "cm"),
          title_gp = gpar(fontsize = 10),
          labels_gp = gpar(fontsize = 8),
          at = c(-2, 0, 2),
          padding = unit(c(2, 2, 2, 2), "mm")
        ),
        row_title = "Proteins",
        column_title = "Samples by Sex & Treatment",
        row_title_gp = gpar(fontsize = 8),
        column_title_gp = gpar(fontsize = 8),
        use_raster = FALSE  # Explicitly set to disable rasterization warning
)
