# Enterococcus-phage-levodopa analysis

Reproducibility package for the main figures of the study **“Precision phage cocktail targeting *Enterococcus faecalis* restores levodopa responsiveness in Parkinson’s disease.”**

The repository is organized as one folder per main figure, following the format requested by *iMeta*. Each figure folder contains an R Markdown (`.Rmd`) workflow together with the corresponding analysis scripts and source/processed input data required for the R-generated panels.

## Repository structure

- `fig1/fig1.Rmd` — computational and statistical analyses used for Figure 1.
- `fig2/fig2.Rmd` — computational and statistical analyses used for Figure 2.

## Figure 1

The Figure 1 workflow includes Enterococcus and Enterococcus-associated phage alpha-diversity analyses, differential abundance analyses, clinical correlation analyses, cross-omics correlation analysis, and LOOCV ROC modeling. Figure 1E uses mutually exclusive significance tiers: `p < 0.05 and q >= 0.10`, `p < 0.05 and 0.05 <= q < 0.10`, and `p < 0.05 and q < 0.05`.

Clinical metadata included in this repository are restricted to the de-identified variables required to reproduce the reported main-figure analyses.

## Figure 2

The Figure 2 workflow includes striatal metabolomics, fecal 16S rRNA community analyses, differential species analysis, and integrated bacteria–metabolite correlation analysis. For the metabolite differential-analysis workflow used in the figure, non-detects/zero values are replaced by half of the minimum positive value for each metabolite, followed by log2 transformation,Wilcoxon rank-sum testing, pooled-variance Fisher's LSD pairwise comparisons, and Benjamini–Hochberg adjustment within each pairwise contrast.

Figure 2M evaluates Spearman correlations across the 15 paired samples from the GR, MR, and MR + Phage groups, with BH correction across the 160 bacteria–metabolite pairs.

## Running the analyses

Open `fig1/fig1.Rmd` or `fig2/fig2.Rmd` in RStudio and knit the document, or execute the corresponding `.R` scripts individually. All analysis paths in the supplied scripts are relative to the corresponding figure directory.

Analyses were performed primarily in R 4.5.2. Required packages are loaded within the scripts and include, depending on the panel, `data.table`, `tidyverse`, `readxl`, `vegan`, `phyloseq`, `ggplot2`, `pheatmap`, `pROC`, `logistf`, `igraph`, `ggraph`, and related packages.

## Public data repositories

Raw sequencing and metabolomics datasets are deposited in the repositories cited in the associated manuscript. The GitHub repository contains the figure-level data and scripts needed for the main-figure reproducibility package.

## License

Code in this repository is released under the MIT License.
