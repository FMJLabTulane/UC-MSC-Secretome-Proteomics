# UC-MSC-Secretome-Proteomics
Secretome Proteomics analysis for male and female UC-MSC donors that were subjected to different sex steroids and inhibitors.
# 🧬 MSC Secretome Proteomics Analysis

This project investigates how sex steroids and their inhibitors affect the secretome of male and female human umbilical cord-derived mesenchymal stem cells (UC-MSCs). Our goal is to uncover key secreted proteins and pathways associated with cellular senescence, metabolic regulation, and regenerative potential.

---

## 📁 Project Structure

```
MSC-Secretome-Analysis/
├── data/                     # Raw and processed iBAQ data from MS
├── scripts/                  # R scripts for analysis and visualization
├── results/                  # Output files: plots, enrichment tables, PCA, etc.
├── README.md                 # Project overview
└── requirements.txt          # R package dependencies
```

---

## 🧪 Experimental Design

* **Cell Type**: Human UC-MSCs (male and female donors)
* **Treatments**:

  * Estradiol (E2)
  * Dihydrotestosterone (DHT)
  * Testosterone (T)
  * Testosterone + Anastrozole
  * Testosterone + 5α-reductase inhibitors (Dutasteride, Finasteride)
  * Vehicle control (EtOH)
* **Secretome Collection**: Serum-free conditioned media collected after 5–7 days of treatment
* **Proteomics**: Mass spectrometry (LC-MS/MS) and iBAQ quantification

---

## 🧾 Data Analysis Pipeline

Analysis was done using R, based on iBAQ values:

1. **Data Preprocessing**:

   * Cleaning, filtering low-abundance proteins
   * Log2 transformation and normalization

2. **PCA**:

   * Identify sample outliers and treatment clustering

3. **Differential Expression**:

   * Linear modeling with `limma`
   * Pairwise contrasts between treatment and control

4. **Functional Enrichment**:

   * Gene Ontology (GO: BP, CC, MF) and KEGG
   * Conducted using `gprofiler2`

5. **Visualization**:

   * Volcano plots
   * Heatmaps of significant DEPs
   * Enrichment plots

---

## 📊 Key Findings

* **Sex-Specific Secretome Profiles** under steroid stimulation
* **Steroid- and inhibitor-dependent alterations** in ECM organization, adhesion, metabolic regulation
* Identification of **SASP-related** protein signatures under DHT and E2 stimulation

---

## 📦 Dependencies

Install required R packages:

```r
install.packages(c("tidyverse", "limma", "gprofiler2", "pheatmap", "ggplot2", "FactoMineR", "factoextra"))
```

---

## 🔁 Reproducibility

To reproduce the full analysis:

```bash
Rscript scripts/run_analysis_pipeline.R
```

*Make sure your `data/` folder contains properly named iBAQ files.*

---

## 📚 Citation / Acknowledgment

If you use this code or workflow, please cite:

```
Nazmul et al., MSC Secretome and Sex Steroid Influence on Senescence, 2025 (preprint/in preparation)
```

Analysis scripts assisted by OpenAI ChatGPT (May 2025, GPT-4.5).

---

## 📬 Contact

**Nazmul**, Postdoctoral Fellow
Mauvais-Jarvis Lab, Tulane University
📧 [your.email@domain.com](mailto:your.email@domain.com)


