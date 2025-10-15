

This repository contains all the materials used for our study **Passive air sampling detects environmental DNA transfer from water into air**. Here, you’ll find everything needed to reproduce analyses, review results, and explore the project — including code, data, model outputs, figures, and the final manuscript.

| Folder | Description |
|:-------|:-------------|
| **Code/** | Contains all the R and Stan scripts used to load, transform, and analyze the data. These scripts include preprocessing steps, the Stan joint model, and plotting routines for all figures. |
| **Data/** | Raw, unprocessed data collected during the field and laboratory work. No transformations applied here — this is the original dataset used as model input. |
| **Output/** | Contains all the .RDS outputs (one in binary form) from the Stan joint model (posterior distributions). |
| **Plots/** | All figures, visualizations, and graphical summaries produced from the R scripts (and beyond). These correspond to manuscript figures and supplementary plots. |
| **Manuscript/** | The final (and older version in Recursive), formatted version of the manuscript submitted for publication (LaTeX and Word). |
| **Review/** | Includes peer-review comments. |
| **Other/** | Mischelleanous documents |

1. **Open the R project**  
   Launch `Fish_can_fly.RProj` in RStudio to access the environment.
2. **Run the scripts**  
   - Data cleaning and transformation scripts are in `Code/`.  
   - The Stan model can be executed directly from the provided R scripts.  
   - Figures are generated automatically after the model runs.  
3. **Explore outputs**  
   Model results and posterior distributions are stored in `Output/`, while generated figures can be found in `Plots/`.
