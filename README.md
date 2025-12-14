# Supplementary Figures: Relevance of DNA Tridimensional Shape in RNA:DNA:DNA Triple Helix Formation

This repository contains the R code required to reproduce the supplementary figures (Supplementary Figures 4 and 5) presented in the manuscript. The analysis correlates the predicted triplex stability (3plex score) with DNA 3D shape features (Helix Twist, Minor Groove Width, Propeller Twist, and Roll).

## 📂 Repository Structure
Ensure your directory contains the following folders with the pre-calculated matrix files before running the script:
* `Random_Negatives/`
* `cCRE_Balanced/`
* `Biosample_Specific/`

## 🛠️ Prerequisites
* **R** (version 4.0.0 or higher recommended)
* The script automatically checks for and installs the following required R packages from CRAN if they are missing:
    * `tidyverse`
    * `ggpubr`
    * `data.table`
    * `hexbin`
    * `ggpointdensity`

## 🚀 How to Run
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/molinerisLab/lncRNA_3plex_DNAshape.git
    ```
2.  **Navigate to the folder:**
    ```bash
    cd lncRNA_3plex_DNAshape
    ```
3.  **Run the analysis script:**
    You can run the script directly from the terminal:
    ```bash
    Rscript script.R
    ```
    *Alternatively, you can open `script.R` in RStudio and click "Source".*

## 📊 Output Description
The script processes three different control datasets (`Random_Negatives`, `cCRE_Balanced`, and `Biosample_Specific`) and generates the following files for each:

### **1. Supplementary Figure 5: Global Density Correlation**
* **Filename:** `RegionPlot_OVERALL_DENSITY.pdf` (and `.png`)
* **Description:** A density-colored scatterplot visualizing the relationship between DNA shape features (x-axis) and triplex stability scores (y-axis) across all genomic binding sites.
* **Visualization:** Points are colored by local density (Yellow = high density, Purple = low density) to visualize the core distribution of stable triplexes.

### **2. Supplementary Figure 4: Faceted lncRNA Analysis**
* **Filename:** `RegionPlot_FACETED_Top10.pdf` (and `.png`)
* **Description:** A detailed breakdown of the correlation for each individual lncRNA.
* **Visualization:** * **Yellow points:** Represent the **top 10%** highest-stability binding sites for that specific lncRNA.
    * **Gray points:** Represent the remaining 90% (lower stability) sites.
    * **Trend lines:** Linear regression fits for each lncRNA.

### **3. Statistical Tables**
* **Filename:** `RegionPlot_FACETED_cor_table_Top10.csv`
* **Description:** A CSV file containing the Pearson correlation coefficients (R) and p-values for every lncRNA-shape pair.

---
