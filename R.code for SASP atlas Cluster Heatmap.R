# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

############### Female DHT DOWN #####################

# 1. Read your data
# (Update your file path if needed)
Female_DHT_Down <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_DHT/DOWN/SASP_Atlas_Female DHT_Down.csv")

# 2. Reshape data into long format
long_Female_DHT_Down <- Female_DHT_Down %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_DHT_Down <- long_Female_DHT_Down %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_DHT_Down <- heatmap_Female_DHT_Down %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_DHT_Down <- heatmap_Female_DHT_Down %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_DHT_Down[is.na(heatmap_Female_DHT_Down)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Female_DHT_Down, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_DHT/DOWN/Filtered_SASP_Atlas_Female_DHT_Down.csv")

# Calculate your global range
min_value <- min(heatmap_Female_DHT_Down, na.rm = TRUE)
max_value <- max(heatmap_Female_DHT_Down, na.rm = TRUE)

# Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)

# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_DHT_Down),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 2,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female DHT DOWN_Overlapping SASP Proteins",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_DHT_Down, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_DHT/DOWN/Filtered_SASP_Atlas_Female_DHT_DOWN.csv")


##############  Female DHT UP #####
# 1. Read your data
# (Update your file path if needed)
Female_DHT_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_DHT/UP/SASP_Atlas_Female DHT_UP.csv")

# 2. Reshape data into long format
long_Female_DHT_UP <- Female_DHT_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_DHT_UP <- long_Female_DHT_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_DHT_UP <- heatmap_Female_DHT_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_DHT_UP <- heatmap_Female_DHT_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_DHT_UP[is.na(heatmap_Female_DHT_UP)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Female_DHT_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_DHT/UP/Filtered_SASP_Atlas_Female_DHT_UP.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_DHT_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 4,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female DHT UP_Overlapping SASP Proteins",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_DHT_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_DHT/UP/Filtered_SASP_Atlas_Female_DHT_UP.csv")

##############  Female E2 DOWN #####
# 1. Read your data
# (Update your file path if needed)
Female_E2_DOWN <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_E2/DOWN/SASP_Atlas_Female E2_Down.csv")

# 2. Reshape data into long format
long_Female_E2_DOWN <- Female_E2_DOWN %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_E2_DOWN <- long_Female_E2_DOWN %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_E2_DOWN <- heatmap_Female_E2_DOWN %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_E2_DOWN <- heatmap_Female_E2_DOWN %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_E2_DOWN[is.na(heatmap_Female_E2_DOWN)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Female_E2_DOWN, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_E2/DOWN/Filtered_SASP_Atlas_Female_E2_DOWN.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_E2_DOWN),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 4,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female E2 DOWN_Overlapping SASP Proteins",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_E2_DOWN, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_E2/DOWN/Filtered_SASP_Atlas_Female_E2_DOWN.csv")

##############  Female E2 UP #####
# 1. Read your data
# (Update your file path if needed)
Female_E2_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_E2/UP/SASP_Atlas_Female E2 UP.csv")

# 2. Reshape data into long format
long_Female_E2_UP <- Female_E2_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_E2_UP <- long_Female_E2_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_E2_UP <- heatmap_Female_E2_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_E2_UP <- heatmap_Female_E2_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_E2_UP[is.na(heatmap_Female_E2_UP)] <- 0


# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_E2_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 2,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female E2 UP_Overlapping SASP Proteins",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_E2_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_E2/UP/Filtered_SASP_Atlas_Female_E2_UP.csv")

##############  Female T DOWN #####
# 1. Read your data
# (Update your file path if needed)
Female_T_DOWN <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_T/DOWN/SASP_Atlas_Female T DOWN.csv")

# 2. Reshape data into long format
long_Female_T_DOWN <- Female_T_DOWN %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_T_DOWN <- long_Female_T_DOWN %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_T_DOWN <- heatmap_Female_T_DOWN %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_T_DOWN <- heatmap_Female_T_DOWN %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_T_DOWN[is.na(heatmap_Female_T_DOWN)] <- 0


# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_T_DOWN),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female T DOWN_Overlapping SASP Proteins",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_T_DOWN, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_T/DOWN/Filtered_SASP_Atlas_Female_T_DOWN.csv")

##############  Female T UP #####
# 1. Read your data
# (Update your file path if needed)
Female_T_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_T/UP/SASP_Atlas_Female T UP.csv")

# 2. Reshape data into long format
long_Female_T_UP <- Female_T_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_T_UP <- long_Female_T_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_T_UP <- heatmap_Female_T_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_T_UP <- heatmap_Female_T_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_T_UP[is.na(heatmap_Female_T_UP)] <- 0


# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_T_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female T UP_Overlapping SASP Proteins",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_T_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_T/UP/Filtered_SASP_Atlas_Female_T_UP.csv")

##############  Female TA UP #####
# 1. Read your data
# (Update your file path if needed)
Female_TA_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_TA/UP/SASP_Atlas_Female TA UP.csv")

# 2. Reshape data into long format
long_Female_TA_UP <- Female_TA_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_TA_UP <- long_Female_TA_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_TA_UP <- heatmap_Female_TA_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_TA_UP <- heatmap_Female_TA_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_TA_UP[is.na(heatmap_Female_TA_UP)] <- 0


# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_TA_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female T+Ana UP_Overlapping SASP Proteins",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_TA_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_TA/UP/Filtered_SASP_Atlas_Female_TA_UP.csv")

##############  Female TA DOWN #####
# 1. Read your data
# (Update your file path if needed)
Female_TA_DOWN <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_TA/DOWN/SASP_Atlas_Female TA DOWN.csv")

# 2. Reshape data into long format
long_Female_TA_DOWN <- Female_TA_DOWN %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_TA_DOWN <- long_Female_TA_DOWN %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_TA_DOWN <- heatmap_Female_TA_DOWN %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_TA_DOWN <- heatmap_Female_TA_DOWN %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_TA_DOWN[is.na(heatmap_Female_TA_DOWN)] <- 0


# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_TA_DOWN),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female T+Ana DOWN_Overlapping SASP Proteins",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_TA_DOWN, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_TA/DOWN/Filtered_SASP_Atlas_Female_TA_DOWN.csv")

##############  Female TS UP #####
# 1. Read your data
# (Update your file path if needed)
Female_TS_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_TS/UP/SASP_Atlas_Female TS UP.csv")

# 2. Reshape data into long format
long_Female_TS_UP <- Female_TS_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_TS_UP <- long_Female_TS_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_TS_UP <- heatmap_Female_TS_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_TS_UP <- heatmap_Female_TS_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_TS_UP[is.na(heatmap_Female_TS_UP)] <- 0


# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_TS_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female T+5α-RIs UP_Overlapping SASP Proteins",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_TS_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_TS/UP/Filtered_SASP_Atlas_Female_TS_UP.csv")

##############  Female TS DOWN #####
# 1. Read your data
# (Update your file path if needed)
Female_TS_DOWN <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_TS/DOWN/SASP_Atlas_Female TS Down.csv")

# 2. Reshape data into long format
long_Female_TS_DOWN <- Female_TS_DOWN %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Female_TS_DOWN <- long_Female_TS_DOWN %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Female_TS_DOWN <- heatmap_Female_TS_DOWN %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Female_TS_DOWN <- heatmap_Female_TS_DOWN %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Female_TS_DOWN[is.na(heatmap_Female_TS_DOWN)] <- 0


# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Female_TS_DOWN),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Female T+5α-RI DOWN_Overlapping SASP Protein",
  angle_col = 45
)

# 8. Save the filtered data
write.csv(heatmap_Female_TS_DOWN, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Female_TS/DOWN/Filtered_SASP_Atlas_Female_TS_DOWN.csv")



############### Male DHT DOWN #####################

# 1. Read your data
# (Update your file path if needed)
Male_DHT_Down <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_DHT/DOWN/SASP_Atlas_Male DHT Down.csv")

# 2. Reshape data into long format
long_Male_DHT_Down <- Male_DHT_Down %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_DHT_Down <- long_Male_DHT_Down %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_DHT_Down <- heatmap_Male_DHT_Down %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_DHT_Down <- heatmap_Male_DHT_Down %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_DHT_Down[is.na(heatmap_Male_DHT_Down)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_DHT_Down, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_DHT/DOWN/Filtered_SASP_Atlas_Male_DHT_Down.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_DHT_Down),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 2,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male DHT DOWN_Overlapping SASP Proteins",
  angle_col = 45
)



##############  Male DHT UP #####
# 1. Read your data
# (Update your file path if needed)
Male_DHT_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_DHT/UP/SASP Atlas_Male DHT UP.csv")

# 2. Reshape data into long format
long_Male_DHT_UP <- Male_DHT_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_DHT_UP <- long_Male_DHT_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_DHT_UP <- heatmap_Male_DHT_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_DHT_UP <- heatmap_Male_DHT_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_DHT_UP[is.na(heatmap_Male_DHT_UP)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_DHT_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_DHT/UP/Filtered_SASP_Atlas_Male_DHT_UP.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down),
  as.matrix(heatmap_Male_DHT_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_DHT_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 2,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male DHT UP_Overlapping SASP Proteins",
  angle_col = 45
)

##############  Male E2 DOWN #####
# 1. Read your data
# (Update your file path if needed)
Male_E2_Down <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_E2/DOWN/SASP_Atlas_Male E2 Down.csv")

# 2. Reshape data into long format
long_Male_E2_Down <- Male_E2_Down %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_E2_Down <- long_Male_E2_Down %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_E2_Down <- heatmap_Male_E2_Down %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_E2_Down <- heatmap_Male_E2_Down %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_E2_Down[is.na(heatmap_Male_E2_Down)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_E2_Down, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_E2/DOWN/Filtered_SASP_Atlas_Male_E2_Down.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down),
  as.matrix(heatmap_Male_DHT_UP),
  as.matrix(heatmap_Male_E2_Down)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_E2_Down),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 2,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male E2 Down_Overlapping SASP Proteins",
  angle_col = 45
)

##############  Male E2 UP #####
# 1. Read your data
# (Update your file path if needed)
Male_E2_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_E2/UP/SASP_Atlas_Male E2 UP.csv")

# 2. Reshape data into long format
long_Male_E2_UP <- Male_E2_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_E2_UP <- long_Male_E2_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_E2_UP <- heatmap_Male_E2_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_E2_UP <- heatmap_Male_E2_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_E2_UP[is.na(heatmap_Male_E2_UP)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_E2_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_E2/UP/Filtered_SASP_Atlas_Male_E2_UP.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down),
  as.matrix(heatmap_Male_DHT_UP),
  as.matrix(heatmap_Male_E2_Down),
  as.matrix(heatmap_Male_E2_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_E2_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 2,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male E2 UP_Overlapping SASP Proteins",
  angle_col = 45
)

##############  Male T DOWN #####
# 1. Read your data
# (Update your file path if needed)
Male_T_Down <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_T/DOWN/SASP_Atlas_Male T Down.csv")

# 2. Reshape data into long format
long_Male_T_Down <- Male_T_Down %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_T_Down <- long_Male_T_Down %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_T_Down <- heatmap_Male_T_Down %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_T_Down <- heatmap_Male_T_Down %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_T_Down[is.na(heatmap_Male_T_Down)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_T_Down, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_T/DOWN/Filtered_SASP_Atlas_Male_T_Down.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down),
  as.matrix(heatmap_Male_DHT_UP),
  as.matrix(heatmap_Male_E2_Down),
  as.matrix(heatmap_Male_E2_UP),
  as.matrix(heatmap_Male_T_Down)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_T_Down),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male T Down_Overlapping SASP Proteins",
  angle_col = 45
)
##############  Male T UP #####
# 1. Read your data
# (Update your file path if needed)
Male_T_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_T/UP/SASP_Atlas_Male T UP.csv")

# 2. Reshape data into long format
long_Male_T_UP <- Male_T_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_T_UP <- long_Male_T_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_T_UP <- heatmap_Male_T_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_T_UP <- heatmap_Male_T_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_T_UP[is.na(heatmap_Male_T_UP)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_T_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_T/UP/Filtered_SASP_Atlas_Male_T_UP.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down),
  as.matrix(heatmap_Male_DHT_UP),
  as.matrix(heatmap_Male_E2_Down),
  as.matrix(heatmap_Male_E2_UP),
  as.matrix(heatmap_Male_T_Down),
  as.matrix(heatmap_Male_T_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_T_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male T UP_Overlapping SASP Proteins",
  angle_col = 45
)
##############  Male TA UP #####
# 1. Read your data
# (Update your file path if needed)
Male_TA_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_TA/UP/SASP_Atlas_UP Male TA.csv")

# 2. Reshape data into long format
long_Male_TA_UP <- Male_TA_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_TA_UP <- long_Male_TA_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_TA_UP <- heatmap_Male_TA_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_TA_UP <- heatmap_Male_TA_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_TA_UP[is.na(heatmap_Male_TA_UP)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_TA_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_TA/UP/Filtered_SASP_Atlas_Male_TA_UP.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down),
  as.matrix(heatmap_Male_DHT_UP),
  as.matrix(heatmap_Male_E2_Down),
  as.matrix(heatmap_Male_E2_UP),
  as.matrix(heatmap_Male_T_Down),
  as.matrix(heatmap_Male_T_UP),
  as.matrix(heatmap_Male_T_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_TA_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 2,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male T+Ana UP_Overlapping SASP Proteins",
  angle_col = 45
)
##############  Male TA DOWN #####
# 1. Read your data
# (Update your file path if needed)
Male_TA_Down <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_TA/DOWN/SASP_Atlas_Male TA Down.csv")

# 2. Reshape data into long format
long_Male_TA_Down <- Male_TA_Down %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_TA_Down <- long_Male_TA_Down %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_TA_Down <- heatmap_Male_TA_Down %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_TA_Down <- heatmap_Male_TA_Down %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_TA_Down[is.na(heatmap_Male_TA_Down)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_TA_Down, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_TA/DOWN/Filtered_SASP_Atlas_Male_TA_DOWN.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down),
  as.matrix(heatmap_Male_DHT_UP),
  as.matrix(heatmap_Male_E2_Down),
  as.matrix(heatmap_Male_E2_UP),
  as.matrix(heatmap_Male_T_Down),
  as.matrix(heatmap_Male_T_UP),
  as.matrix(heatmap_Male_T_UP),
  as.matrix(heatmap_Male_TA_Down)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_TA_Down),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 2,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male T+Ana Down_Overlapping SASP Proteins",
  angle_col = 45
)
##############  Male TS DOWN #####
# 1. Read your data
# (Update your file path if needed)
Male_TS_DOWN <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_TS/DOWN/SASP_Atlas_Male TS_DOWN.csv")

# 2. Reshape data into long format
long_Male_TS_DOWN <- Male_TS_DOWN %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_TS_DOWN <- long_Male_TS_DOWN %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_TS_DOWN <- heatmap_Male_TS_DOWN %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_TS_DOWN <- heatmap_Male_TS_DOWN %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_TS_DOWN[is.na(heatmap_Male_TS_DOWN)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_TS_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_TS/DOWN/Filtered_SASP_Atlas_Male_TS_DOWN.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down),
  as.matrix(heatmap_Male_DHT_UP),
  as.matrix(heatmap_Male_E2_Down),
  as.matrix(heatmap_Male_E2_UP),
  as.matrix(heatmap_Male_T_Down),
  as.matrix(heatmap_Male_T_UP),
  as.matrix(heatmap_Male_T_UP),
  as.matrix(heatmap_Male_TA_Down),
  as.matrix(heatmap_Male_TS_DOWN)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_TS_DOWN),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male T+5αRI DOWN_Overlapping SASP Proteins",
  angle_col = 45
)
##############  Male TS UP #####
# 1. Read your data
# (Update your file path if needed)
Male_TS_UP <- read_csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_TS/UP/SASP_Atlas_Male TS_UP.csv")

# 2. Reshape data into long format
long_Male_TS_UP <- Male_TS_UP %>%
  pivot_longer(
    cols = c(`Fibro IR`, `Fibro RAS`, `Fibro ATV`, `Epi IR`, `Exosome IR`, `Exosome RAS`),
    names_to = "SASP_Category",
    values_to = "Expression"
  )

# 3. Prepare data matrix for heatmap
heatmap_Male_TS_UP <- long_Male_TS_UP %>%
  pivot_wider(names_from = SASP_Category, values_from = Expression) %>%
  column_to_rownames(var = "Genes")

# 4. Remove the 'Epi IR' column
heatmap_Male_TS_UP <- heatmap_Male_TS_UP %>%
  select(-`Epi IR`)

# 5. Remove proteins with all NA

heatmap_Male_TS_UP <- heatmap_Male_TS_UP %>%
  filter(!(if_all(everything(), is.na)))

# 6. Replace NA with 0 or you can choose to replace with row mean if you prefer
heatmap_Male_TS_UP[is.na(heatmap_Male_TS_UP)] <- 0

# 7. Save the filtered data
write.csv(heatmap_Male_TS_UP, "C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/UC-MSC project/Raw data&analysis/Secretome-Proteomics/Solo Analysis by NH/SASP atlas analysis/Male_TS/UP/Filtered_SASP_Atlas_Male_TS_UP.csv")

# 1. Combine all datasets into one temporary matrix
combined_data <- rbind(
  as.matrix(heatmap_Female_DHT_Down),
  as.matrix(heatmap_Female_DHT_UP),
  as.matrix(heatmap_Female_E2_DOWN),
  as.matrix(heatmap_Female_E2_UP),
  as.matrix(heatmap_Female_T_DOWN),
  as.matrix(heatmap_Female_T_UP),
  as.matrix(heatmap_Female_TA_UP),
  as.matrix(heatmap_Female_TA_DOWN),
  as.matrix(heatmap_Female_TS_UP),
  as.matrix(heatmap_Female_TS_DOWN),
  as.matrix(heatmap_Male_DHT_Down),
  as.matrix(heatmap_Male_DHT_UP),
  as.matrix(heatmap_Male_E2_Down),
  as.matrix(heatmap_Male_E2_UP),
  as.matrix(heatmap_Male_T_Down),
  as.matrix(heatmap_Male_T_UP),
  as.matrix(heatmap_Male_T_UP),
  as.matrix(heatmap_Male_TA_Down),
  as.matrix(heatmap_Male_TS_UP),
  as.matrix(heatmap_Male_TS_UP)
)

# 2. Calculate global min and max from combined data
min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)


# 3. Create consistent breaks
my_breaks <- seq(min_value, max_value, length.out = 100)


# Now plot the heatmap with fixed scale
pheatmap(
  as.matrix(heatmap_Male_TS_UP),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 4,
  fontsize_col = 10,
  color = colorRampPalette(c("lightblue", "blue", "red"))(length(my_breaks) - 1),
  breaks = my_breaks,
  main = "Male T+5αRI UP_Overlapping SASP Proteins",
  angle_col = 45
)
