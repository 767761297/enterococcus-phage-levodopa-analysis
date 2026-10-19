# Enterococcus-phage-levodopa analysis

Reproducibility package for the main figures of the study **“Precision phage cocktail restores levodopa response in Parkinson’s disease.”**

The repository is organized with figure-specific analysis/data folders and a separate `scripts/` directory for standalone R scripts, following the organization requested during *iMeta* production.

## Repository structure

- `fig1/fig1.Rmd` — executable R Markdown workflow and source/processed data for the R-generated panels of Figure 1.
- `fig2/fig2.Rmd` — executable R Markdown workflow and source/processed data for the R-generated panels of Figure 2.
- `scripts/fig1/` — standalone R scripts corresponding to the Figure 1 analyses.
- `scripts/fig2/` — standalone R scripts corresponding to the Figure 2 analyses.
- `scripts/README.md` — script inventory and execution instructions.

## Figure 1

The Figure 1 workflow includes Enterococcus and Enterococcus-associated phage alpha-diversity analyses, differential abundance analyses, clinical correlation analyses, cross-omics correlation analysis, and LOOCV ROC modeling.

Figure 1E uses mutually exclusive significance tiers: `p < 0.05 and q >= 0.10`, `p < 0.05 and 0.05 <= q < 0.10`, and `p < 0.05 and q < 0.05`.

Clinical metadata included in this repository are restricted to the de-identified variables required to reproduce the reported main-figure analyses.

## Figure 2

The Figure 2 workflow includes striatal metabolomics, fecal 16S rRNA community analyses, differential species analysis, and integrated bacteria–metabolite correlation analysis.

For the metabolite differential-analysis workflow used in the final figure, non-detects/zero values are replaced by half of the minimum positive value for each metabolite, followed by log2 transformation, one-way ANOVA, pooled-variance Fisher's LSD pairwise comparisons, and Benjamini–Hochberg adjustment within each pairwise contrast.

Figure 2M evaluates Spearman correlations across the 15 paired samples from the GR, MR, and MR + Phage groups, with BH correction across the 160 bacteria–metabolite pairs.

## Running the analyses

For the figure-level reproducible workflows, open `fig1/fig1.Rmd` or `fig2/fig2.Rmd` in RStudio and knit from the corresponding figure directory.

The standalone R scripts are stored separately under `scripts/`. Because their input/output paths are relative to the corresponding figure directory, run them with that figure directory as the working directory. For example:

```bash
cd fig1
Rscript ../scripts/fig1/Fig1B-1.R
```

or

```bash
cd fig2
Rscript ../scripts/fig2/Fig2G.R
```

See `scripts/README.md` for the complete script inventory and Figure 2 execution note.

Analyses were performed primarily in R 4.5.2. Required packages are loaded within the scripts and include, depending on the panel, `data.table`, `tidyverse`, `readxl`, `vegan`, `phyloseq`, `ggplot2`, `pheatmap`, `pROC`, `logistf`, `igraph`, `ggraph`, and related packages.

## Public data repositories

Raw sequencing and metabolomics datasets are deposited in the public repositories cited in the associated manuscript. This GitHub repository contains the figure-level data and scripts needed for the main-figure reproducibility package.

## License

Code in this repository is released under the MIT License.
