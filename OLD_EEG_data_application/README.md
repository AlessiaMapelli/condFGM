# EEG Data Application for condFGM

This folder contains the real-world application of the conditional Functional Graphical Model (condFGM) method to electroencephalography (EEG) data for studying genetic predisposition to alcoholism.

## Dataset Description

The EEG data analyzed in this application was obtained from the repository associated with the work by Zhao, B., Wang, Y. S. & Kolar, M. (2022), "FuDGE: A method to estimate a functional differential graph in a high-dimensional setting", *Journal of Machine Learning Research* 23(82), 1–82.

**Data Source**: https://github.com/boxinz17/FuDGE/tree/master/EEG/

### Data Files

This application uses two versions of the EEG dataset:

#### 1. `alco_array.Rdata` (Raw Files)
- Contains the original, unprocessed EEG recordings
- Includes raw signal data from all 64 electrodes
- Suitable for custom preprocessing pipelines
- Format: R data array

#### 2. `alco_filtered_array.Rdata` (Preprocessed Data)
- Contains filtered and preprocessed EEG signals
- Ready for direct analysis with condFGM
- Preprocessing likely includes noise reduction and artifact removal
- Format: R data array

## Application Overview

The application to real word data follows the computational template instructions available in the README. 
Specifically:

###### 1. Prepare Your Data
- **`scores_df`**:
- - **`covariates_df`**:

###### 2. Configure Parameters
```yaml
# Essential Parameters
input_path: "your_data.RData"
output_path: "results/"
name_output: "my_analysis"
n_nodes: 10                    # Number of brain nodes
n_basis_for_dim_reduction: 5   # Basis functions per node

# Model Parameters
L: 100                         # Number of lambda values to be tested
K: 5                          # K-fold cross-validation
thres_ctrl: [0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]  # Threshold grid
p_rand_lam: 0.5               # Proportion of lambda values to test
p_rand_thr: 1                 # Proportion of thresholds to test
type: OR                      # Edge symmetrization method. Options: OR (Edge exists if detected in either direction), AND (Edge exists only if detected in both directions)

verbose: TRUE                 #Enable progress messages. Options: TRUE (Recommended for monitoring), FALSE
```
###### 3. Execute Sequential Processing
```bash
# Run complete analysis on single machine
Rscript Script_sequential.R config_file.yaml
```

###### 4. Process Results

```bash
# Generate final adjacency matrices and visualizations
Rscript Results_evaluation.R config_file.yaml
```

**Note**: This application demonstrates the practical utility of condFGM for neuroscience research and provides a benchmark for comparing functional differential graph estimation methods on real neurological data.

## Output Files

###### 5. Plot final results
Run Eval_FudgeDataset_topomap.R to visualize the differential network.

