# LANSA: Laminar Anatomy and Neural Systems Atlas

Analysis pipeline for:

> Mulvey AG, Lichtenfeld MJ, Nejat H, Xiong Y, Carlson BM, Mitchell BA, Mendoza-Halliday D, Westerberg JA, Maier A, Duque A, Bastos AM. **Conserved Laminar-Cellular Architecture Predicts the Spectrolaminar Motif.**

## Contents

| File | Description |
|---|---|
| `LANSA_Analysis_Pipeline.m` | Complete analysis pipeline: cross-area conservation, PCA, Procrustes analysis, PLSR and the schizophrenia perturbation (Figures 2B, 3–6, S3–S5, S7) |
| `LANSA_processed_cyto.mat` | Preprocessed cytoarchitecture: laminar density of PV, CB, CR and NG cells in 11 areas, in the standardized laminar space (21 bins), plus the laminar boundaries of each area |
| `functions/` | `buildFeatureMatrix` (PLSR predictors), `permTestMeanR` (conservation permutation test), `cellDensityPlotEditor` (interactive density editor), and plotting helpers |

## Spectrolaminar data

The spectrolaminar (LFP) data are **not included**. They must be added from Mendoza-Halliday et al. (2024), *Nat Neurosci* 27, 547–560 (https://doi.org/10.1038/s41593-023-01554-7).

Save them in the repository folder as `spectrolaminar_data.mat`, containing the variable `ephys_data`. This is relative power with dimensions 21 laminar bins × 150 frequencies (1–150 Hz) × 11 areas. It must use the standardized laminar space (`standard_edges`) and the area order (`included_areas`) of `LAP_pub_data.mat`. The pipeline stops with an error if this file is missing.

## Running

Requires MATLAB R2021a or later and the Statistics and Machine Learning Toolbox.

1. Add `spectrolaminar_data.mat` to the repository folder.
2. Run `LANSA_Analysis_Pipeline.m`.

Section 10 opens an optional interactive editor for the density profiles. Set `RUN_INTERACTIVE_EDITOR = false` at the top of the script to skip it.

Permutation tests and bootstrap resampling use 5,000 iterations with a fixed seed (`rng(42)`).
