# Tissue-Aware Drug Response Prediction

Code and results for two companion papers on drug response prediction from genomic features:

1. **Asiaee et al. (2026)** "Widespread data leakage inflates performance estimates in cancer drug response prediction." *bioRxiv*. doi: [10.1101/2026.02.05.704016](https://doi.org/10.1101/2026.02.05.704016)

2. **Strauch, Azinfar, Pal, Pua, Long, Coombes, Asiaee (2026)** "Tissue-aware elastic net decomposition reveals shared and lineage-specific drug response biomarkers." *In preparation for Bioinformatics (Oxford).*

## Overview

This repository provides:

- **Data Shared Elastic Net (DSEN)**: A tissue-aware regularized regression that decomposes drug response coefficients into shared (pan-cancer) and tissue-specific blocks.
- **Leakage-free cross-validation**: Correct implementation with per-fold feature screening, contrasted with the leaky approach used by the majority of published methods.
- **Method comparison**: DSEN vs standard elastic net vs TG-LASSO across 265 drugs.
- **Biological validation**: Target recovery, COSMIC enrichment, and systematic gene-level analysis of shared vs tissue-specific features.

## Key Results

- DSEN improves MSE for **92.5%** of 265 drugs (mean 4.95% improvement) under leakage-free evaluation.
- TG-LASSO, the only comparable tissue-aware method, performs **worse** than the tissue-agnostic baseline (-13.8% mean MSE).
- Leave-tissue-out validation confirms shared coefficients generalize to unseen tissues (59.3% win rate).
- The strongest biological pattern: MAPK-pathway drugs recover melanoma-lineage markers (MITF, S100B) in skin tissue blocks, linking to concrete leave-skin-out transfer failures.

## Repository Structure

```
methods/           # Reusable R functions
  00_paths.R       # Path configuration
  01_data_utils.R  # Data loading and preprocessing
  02_cv_utils.R    # Cross-validation (correct SE, leakage-free)
  03_dsen_utils.R  # Data Shared Elastic Net implementation
  04_plot_utils.R  # Visualization
  05_parallel_utils.R

experiments/       # Analysis scripts and R Markdown notebooks
results/           # Generated outputs (CSVs, figures)
data/              # Drug target mappings and reference files
```

## Data Requirements

Molecular feature data (CCLE/DepMap) and drug response data (GDSC) must be obtained from their original sources:

- **CCLE/DepMap**: https://depmap.org
- **GDSC**: https://www.cancerrxgene.org

Create a JSON configuration file at `~/Paths/drugresponseprediction.json`:

```json
{
  "paths": {
    "raw": "/path/to/raw/data",
    "clean": "/path/to/clean/data",
    "scratch": "/path/to/scratch"
  }
}
```

## R Dependencies

```r
install.packages(c("here", "glmnet", "caret", "tidyverse", "rjson", "Matrix"))
```

## Citation

If you use this code, please cite:

```bibtex
@article{Asiaee_2026_leakage,
  title={Widespread data leakage inflates performance estimates in cancer drug response prediction},
  author={Asiaee, Amir and Strauch, Jared and Azinfar, Leila and Pal, Samhita and Pua, Heather H. and Long, James P. and Coombes, Kevin R.},
  journal={bioRxiv},
  year={2026},
  doi={10.1101/2026.02.05.704016}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.

