# Structural Damage Detection with Domain Adaptation

This repository contains MATLAB codes for structural health monitoring and damage detection. The project uses Domain Adaptation methods, specifically Transfer Component Analysis (TCA), combined with the Mahalanobis distance to classify experimental data based on a numerical model.

## 📂 Repository Structure

The project is divided into three main modules:

*   **Module 1: Experimental Extraction (`01_Extracao_Experimental.m` and `01_Extracao_Experimental_Atualizado.m`)**
    *   Processes experimental signals, extracting force and accelerometer information.
    *   Extracts natural frequencies using Power Spectral Density (PSD) via Welch's method or the Modal Circle Fit method from complex FRF estimates.
    *   Allows extractions using a single channel or through multichannel analysis with spatial fusion.
    *   Generates and consolidates the physical database in the `Dados_Alvo.mat` (Target Domain) file.

*   **Module 2: Numerical Generation (`02_Geracao_Numerica.m`)**
    *   Generates the Source Domain by simulating theoretical frequencies from a 4-degree-of-freedom (4-DOF) numerical dynamic model.
    *   Simulates structural conditions through stochastic mapping of the Modulus of Elasticity, operational mass variations, and progressive stiffness loss.
    *   Saves the theoretical simulation database in the `Dados_Origem.mat` file.

*   **Module 3: Adaptation and Classification (`03_Adaptacao_Classificacao.m` and `TCA_linear.m`)**
    *   Loads the unified data from the target (`Dados_Alvo.mat`) and source (`Dados_Origem.mat`) domains.
    *   Applies Transfer Component Analysis (TCA) with a linear kernel to align the data in a new latent space and minimize domain discrepancy.
    *   Employs Mahalanobis classification, where the decision boundary learns exclusively from healthy samples (conditions 1 to 3).
    *   Evaluates predictions against damaged samples (conditions 4 to 9) and calculates final performance using the Macro F1 metric.

## 🚀 How to Run

To ensure the scripts work correctly, follow the execution order below:

1.  **Phase 1:** Run one of the Module 1 scripts (`01_Extracao_Experimental.m` or the updated version) to process the files in the `Data` folder and generate the `Dados_Alvo.mat` file.
2.  **Phase 2:** Run `02_Geracao_Numerica.m` to execute the theoretical stochastic model and generate the `Dados_Origem.mat` file.
3.  **Phase 3:** With both files generated in the previous phases, run `03_Adaptacao_Classificacao.m` to start domain alignment and view classification performance.

*** 
*Note: The `Dados_Alvo.mat` and `Dados_Origem.mat` matrices are mandatory requirements to run Phase 3, and the script will return an error if they are not found in the directory.*

## 📊 Dataset

The experimental data required to run the extraction scripts is not included in this repository. You can download the complete dataset from the following link:

👉 **[Download Dataset Here](https://drive.google.com/drive/folders/1ZKao-fyNRudde7wgwp9wjfajUt0-hwDo?usp=sharing)**

**Instructions:**
1. Download the dataset from the link above.
2. Extract the downloaded file.
3. Place the `Data` folder directly into the root directory of this repository. 
4. Proceed with running Phase 1 (`01_Extracao_Experimental.m`).
