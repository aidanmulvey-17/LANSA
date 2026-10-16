# LANSA: Laminar Anatomy and Neural Systems Atlas

Analysis pipeline for:

> Mulvey AG, Lichtenfeld MJ, Nejat H, Xiong Y, Carlson BM, Mitchell BA, Mendoza-Halliday D, Westerberg JA, Maier A, Duque A, Bastos AM. **Conserved Laminar-Cellular Architecture Predicts the Spectrolaminar Motif.**

## Contents

| File | Description |
|---|---|
| `LANSA_Analysis_Pipeline.m` | Complete analysis pipeline: cross-area conservation, PCA, Procrustes analysis, PLSR and the schizophrenia perturbation (Figures 2B, 3–6, S3–S5, S7) |
| `LANSA_processed_cyto.mat` | Preprocessed cytoarchitecture: laminar density of PV, CB, CR and NG cells in 11 areas, in the standardized laminar space (21 bins), plus the laminar boundaries of each area |
| `spectrolaminar_data.mat` | Preprocessed spectrolaminar data in 11 areas, in the standardized laminar space (21 bins), pre-aligned to histology |
| `functions/` | `buildFeatureMatrix` (PLSR predictors), `permTestMeanR` (conservation permutation test), `cellDensityPlotEditor` (interactive density editor), and plotting helpers |

## Running

Requires MATLAB R2021a or later and the Statistics and Machine Learning Toolbox.

Permutation tests and bootstrap resampling use 5,000 iterations with a fixed seed (`rng(42)`).
