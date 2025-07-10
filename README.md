# 🧬 HPA Cell Line Section (v22)

This repository contains the code and resources necessary to reproduce the main findings of the paper:  
**"Systematic transcriptional analysis of human cell lines for gene expression landscape and tumor representation"**

## 🗂 Repository structure
```bash
.
├── data/        # Processed datasets (access via Zenodo)
├── codes/       # Numbered R scripts used in the analysis
├── results/     # Output plots
├── LICENSE
└── README.md
```
- **data/**: Processed data used in the analyses can be accessed via Zenodo:  
  🔗 [https://doi.org/10.5281/zenodo.7874749](https://doi.org/10.5281/zenodo.7874749)
  
- **codes/**: Contains ordered R scripts (e.g., `umap_visualization.R`, `expression_landscape.R`, etc.) for reproducing figures and results presented in the manuscript. The recommended execution order of the scripts is provided in `main.R`.

## ▶️ How to run

1. Clone this repository:

   ```bash
   git clone https://github.com/jha14/HPA_cell_line.git
   cd HPA_cell_line

3. Download the processed data from Zenodo and place it in the data/ directory.
4. Run the scripts sequentially from the codes/ folder:

   ```bash
   source("codes/theme.R")
   source("codes/hpa_functions.R")
   ...

6. Output figures and tables will be saved to the results/ folder as specified.

## 📄 Citation
If you use this code or data in your research, please cite the original paper:

Han Jin, Cheng Zhang, et al. (2023). Systematic transcriptional analysis of human cell lines for gene expression landscape and tumor representation. Nature Communications. [https://doi.org/10.1038/s41467-023-41132-w]

(*Han Jin and Cheng Zhang contributed equally to this work*)
