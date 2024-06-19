# Code repository for the paper "Detecting Alzheimer’s Disease Stages and Frontotemporal Dementia in Time Courses of Resting-State fMRI Data Using a Machine Learning Approach"

You may find the full text of the paper [here](https://link.springer.com/article/10.1007/s10278-024-01101-1).

## Notes:
* The data used in this study were obtained from the ADNI and NIFD databases. The rsfMRIs used in this study were preprocessed and converted to time course data using CONN. Preprocessed data may be shared by the corresponding author of the paper with researchers who have obtained data access licenses from ADNI and NIFD.
* The "experiments" folder contains two subdirectories: "preprocessing_scripts" and "analysis_scripts". Various model architectures (Transformers, LSTM-based NNs, etc.), preprocessing strategies (e.g., implementation of wavelet transformation), and hyperparameter optimization strategies (e.g., use of Weights & Biases) were tested before reaching the final methodologies used in the paper.
* The "final" folder contains the final scripts that were used to preprocess and analyze the data and generate the results and figures that were used in the paper.
