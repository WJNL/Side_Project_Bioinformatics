# Classical Machine Learning and Biologically Informed Feature Representation for Single-Cell Classification

*A comparison of gene-level and pathway-level features on the PBMC 3k dataset*

Code and full explanations: [`Single_cell_ML_project_clean.ipynb`](Single_cell_ML_project_clean.ipynb)

---

## Abstract

Single-cell RNA sequencing (scRNA-seq) gives us gene expression for thousands of individual cells, and machine learning can use it to predict cell types. In this project, I tested two things on the PBMC 3k dataset. First, I compared classical machine learning models (Logistic Regression, Random Forest, SVM and XGBoost) on gene-level features. Second, I replaced the genes with pathway scores built from immune-related Reactome pathways and checked what changes.

The tuned Logistic Regression with 997 genes performed best (test accuracy 0.943, macro F1 0.925). Its main weakness was the separation of CD4 and CD8 T cells. About 16% of the CD8 T cells in the test set were predicted as CD4 T cells by all four models. The pathway-level model with 30 features performed clearly worse (test accuracy 0.778, macro F1 0.697). In simple words: summarising genes into pathway averages made the data easier to read biologically, but it also removed information that the model needed.

**Keywords:** scRNA-seq, cell-type classification, feature selection, Logistic Regression, pathway analysis, Reactome, model interpretability

---

## 1. Introduction

scRNA-seq measures how active each gene is in every single cell. A typical dataset contains thousands of genes per cell, so it is a good case for supervised machine learning. However, gene-level models often use hundreds or thousands of features. They can predict well, but it is difficult to link their decisions to real biological processes.

One idea is to use biological knowledge, for example pathway databases, to build the features. In theory, this could make the model easier to understand. But it is not clear if it keeps the predictive power. This project looks at both sides.

### Research questions

- **RQ1 (gene level):** How effectively can classical machine learning methods classify cell types using gene-level expression features?
- **RQ2 (pathway level):** How does biologically informed feature representation (pathway scores) affect classification performance, and what are the trade-offs compared with gene-level features?

---

## 2. Methods

### 2.1 Dataset

I used the processed PBMC 3k dataset (peripheral blood mononuclear cells) from the `scanpy` package. It has 2,638 cells and 1,838 genes, and eight annotated cell types: CD4 T cells, CD8 T cells, B cells, NK cells, CD14+ monocytes, FCGR3A+ monocytes, dendritic cells and megakaryocytes.

I deliberately selected the **already processed version**. Quality control, normalisation, log transformation, selection of highly variable genes and scaling were already done by the data provider. My focus was on modelling and feature representation, and this choice also saved time. The downside is that I did not build or check the preprocessing myself, and that only highly variable genes are available (this matters later for the pathway coverage).

The classes are strongly imbalanced (Figure 1): there are 1,144 CD4 T cells but only 15 megakaryocytes. For this reason, I used **macro F1** as the main metric. In simple words, macro F1 gives every cell type the same weight, so rare cell types are not ignored like they can be in accuracy.

![Figure 1](images/01_cell_type_distribution.png)

*Figure 1. Number of cells per cell type. The dataset is strongly imbalanced.*

### 2.2 Data split and validation

The data were split into 80% training (2,110 cells) and 20% test (528 cells) with stratified sampling, so that every cell type keeps its proportion. Hyperparameters were tuned with stratified 5-fold cross-validation on the training set only. The test set was used once, for the final evaluation. The same split was used for all models and for the pathway-level analysis.

Evaluation metrics were accuracy, macro F1, weighted F1 and confusion matrices.

### 2.3 RQ1: gene-level classification

- **Feature selection.** An ANOVA F-test (`SelectKBest`) selected the most informative genes. The number of genes (k) was tuned inside the cross-validation pipeline together with the regularisation strength, so the validation folds could not influence the selection. The best result was k = 1,000.
- **Redundancy reduction.** For every pair of selected genes with an absolute correlation above 0.8, one gene was removed. Only three genes were removed, which left **997 genes**.
- **Models.** Logistic Regression, Random Forest, Support Vector Machine (RBF kernel) and XGBoost. Class weights were balanced where possible. XGBoost was used with default settings.

### 2.4 RQ2: pathway-level classification

For the biological representation, I used immune-related pathways from the Reactome database (`Reactome_2022`, loaded with `gseapy`). I selected pathways by keywords (for example immune, interferon, cytokine, antigen, interleukin). For every pathway, I identified the genes that are in the dataset, and I excluded pathways with fewer than five matched genes. This left **30 pathways**.

The pathway score of a cell is the **mean standardised expression** of the matched genes. These 30 scores were the input of a Logistic Regression (with standardisation, balanced class weights and tuned `C`), using the same cells, split and folds as the gene-level model. As a first exploration, I also tested the 50 Hallmark gene sets, but these gave weaker results (test macro F1 0.62), probably because they are not designed for immune cell types.

One important point: these scores summarise the expression of the genes in a pathway. They are **not** direct measurements of pathway activity.

---

## 3. Results

### 3.1 RQ1: gene-level classification

Table 1 shows the test performance of the four models.

**Table 1. Gene-level models (test set, 997 genes)**

| Model | Accuracy | Macro F1 | CD4 recall | CD8 recall |
|---|---:|---:|---:|---:|
| **Logistic Regression (tuned)** | **0.943** | **0.925** | 0.952 | 0.825 |
| SVM (RBF, tuned) | 0.941 | 0.920 | 0.952 | 0.794 |
| XGBoost (default) | 0.930 | 0.893 | 0.952 | 0.778 |
| Random Forest (200 trees) | 0.917 | 0.858 | 0.983 | 0.635 |

The tuned Logistic Regression had the best result on all metrics. Its cross-validated macro F1 was 0.936, and the best regularisation was strong (`C = 0.01`). This makes sense: we have many features and relatively few cells, so a strongly regularised linear model is not easy to overfit. Also, using 1,000 selected genes worked as well as using all 1,838 genes, so many genes did not add much information.

The Random Forest had the highest CD4 recall but the lowest CD8 recall. This suggests that it tends to give uncertain T cells to the large CD4 class.

A note on the SVM: its cross-validated macro F1 (0.823) was much lower than its test macro F1 (0.920). The smallest classes have only two or three cells per validation fold, so the fold scores are unstable. I would not compare these two numbers directly.

### 3.2 Error analysis: CD4 and CD8 T cells

Overall accuracy was high, but the confusion matrix (Figure 2) shows a clear weak point: CD4 and CD8 T cells are mixed up. In the tuned Logistic Regression, 10 CD4 cells were predicted as CD8, and 10 CD8 cells were predicted as CD4.

![Figure 2](images/02_gene_level_confusion_matrix.png)

*Figure 2. Confusion matrix of the tuned gene-level Logistic Regression.*

**Is the problem specific to one algorithm?** No. Out of 63 CD8 T cells in the test set, 42 were never predicted as CD4, but **10 cells (15.9%) were predicted as CD4 by all four models**. So the difficulty is probably in the cells or in the features, not in the algorithm.

**Which genes drive the CD8/CD4 decision?** In a multinomial Logistic Regression, the difference between the CD8 and CD4 coefficients is the weight of each gene in the CD8-versus-CD4 decision. Figure 3 shows the 20 genes with the largest difference. Cytotoxicity-related genes (`CCL5`, `NKG7`, `GZMK`, `CST7`, `GZMA`) push towards CD8, while `TNFRSF4` and `MAL` push towards CD4.

![Figure 3](images/03_cd8_vs_cd4_genes.png)

*Figure 3. Top 20 genes by coefficient difference (CD8 minus CD4). Positive values favour CD8.*

**How do the misclassified cells look?** Figure 4 compares the mean expression of these genes in three groups. The misclassified CD8 cells have low cytotoxic gene expression. For example, the mean `CCL5` value is 1.81 in correctly classified CD8 cells, −0.32 in correctly classified CD4 cells and 0.03 in the misclassified CD8 cells. So, for the model, these cells look more like CD4 cells. The group is not uniform, however: some cells have high `LYAR`, `KLRB1` or `GYG1`.

![Figure 4](images/04_misclassified_cd8_heatmap.png)

*Figure 4. Mean expression of the top 20 discriminative genes in correct CD4 cells, misclassified CD8 cells and correct CD8 cells.*

**A closer look at one cell.** I analysed cell 71, the cell with the strongest preference for CD4 (P(CD4) = 0.874, P(CD8) = 0.046). For each gene, the contribution to the CD8-versus-CD4 score was calculated as expression × (CD8 coefficient − CD4 coefficient). The contributions plus the intercept difference reproduced the model score exactly (absolute error around 10⁻¹⁶), so this explanation is exact for the model. Two patterns appear (Figure 5):

1. Low expression of cytotoxic genes (`CCL5`, `NKG7`, `CTSW`, `CST7`, `GZMA`, `GZMK`) reduces the CD8 score.
2. High expression of genes with CD4-favouring weights (`MAL`, `NT5C3A`, `AQP3`) adds to the push towards CD4.

Some genes support CD8 (for example `LRBA`, +0.333), but not enough. The decision is therefore the sum of many small contributions, not one dominant gene. This analysis explains the *model*; it does not prove that these genes are the biological reason for the error.

![Figure 5](images/05_cell71_contributions.png)

*Figure 5. Gene contributions to the CD8-versus-CD4 score of cell 71 (15 strongest in each direction). Negative values push towards CD4.*

### 3.3 RQ2: pathway-level classification

The Reactome model (30 pathway scores) reached a cross-validated macro F1 of 0.742, a test accuracy of 0.778 and a test macro F1 of 0.697. It recognised clearly distinct cell types quite well (B cells F1 0.93, CD14+ monocytes 0.83, NK cells 0.78). But it was much weaker for similar cell types: CD8 T cells (F1 0.48), FCGR3A+ monocytes (0.52) and dendritic cells (0.32). In the confusion matrix (Figure 6), 31 CD4 T cells were predicted as CD8, 17 CD8 cells as CD4, and 12 more CD8 cells as NK cells. The two monocyte types were also confused.

![Figure 6](images/06_reactome_confusion_matrix.png)

*Figure 6. Confusion matrix of the tuned Reactome pathway-level Logistic Regression.*

### 3.4 Comparison of the two representations

**Table 2. Gene-level vs. pathway-level Logistic Regression**

| Metric | Gene-level (997 features) | Reactome pathway-level (30 features) |
|---|---:|---:|
| Cross-validation macro F1 | 0.936 | 0.742 |
| Test accuracy | 0.943 | 0.778 |
| Test macro F1 | 0.925 | 0.697 |

The gene-level model was better on all metrics.

---

## 4. Discussion

**Classical models for single-cell classification.** Regularised linear models worked best here, and they are also the easiest to interpret. The error analysis was important: an accuracy of 94% hides the fact that one specific pair of cell types (CD4/CD8) causes most of the errors.

**Why did the pathway features perform worse?** A likely reason is information loss. Averaging gives every gene in a pathway the same weight, but genes in one pathway can behave differently. In more technical words, a pathway score is a fixed linear combination of genes, so a Logistic Regression on pathway scores is a restricted version of a Logistic Regression on genes, with equal weights inside each pathway. However, I cannot say that aggregation is the only reason. Other possible factors are:

- **low gene coverage:** the dataset contains only highly variable genes, so on average only about 17% of each pathway's genes are available;
- **keyword-based pathway selection** and **overlapping pathways** (parent and child pathways share many genes);
- the **much lower dimension** (30 vs. 997 features).

I did not test these factors one by one, so their individual effect is unknown.

**Interpretability.** Pathway features are linked to named biological processes, so they have a potential interpretability advantage. But biologically annotated features do not automatically make a model interpretable. The mean score ignores regulation, activation or inhibition, and the different roles of genes in a pathway. Also, I did not test whether pathway-level predictions are easier to explain than gene-level ones. The pathway model is therefore an initial exploration, not a fully interpretable biological model. In my view, the main lesson is that adding biological knowledge needs more than grouping genes by pathway.

---

## 5. Conclusion

1. **RQ1:** Classical machine learning classifies PBMC cell types well. The best model was a tuned Logistic Regression with 997 genes (test accuracy 0.943, macro F1 0.925). The remaining difficulty is the separation of CD4 and CD8 T cells, and the gene-level analysis links it to weaker cytotoxic gene expression in the misclassified CD8 cells.
2. **RQ2:** Reactome pathway scores (30 features) reached a test accuracy of 0.778 and macro F1 of 0.697. This shows a trade-off between keeping detailed gene information and adding biological structure. Simple pathway averaging was not enough to keep the information needed for accurate classification.

---

## 6. Limitations and future work

- **One small dataset, processed version.** Only PBMC 3k was used, without independent validation, and I did not do the raw-data preprocessing myself. Starting from raw counts would be a useful check.
- **No donor information.** If cells from one donor appear in both training and test sets, the results may overestimate generalisation to new donors.
- **Small classes.** The test set contains only 3 megakaryocytes and 7 dendritic cells, so class-level scores and macro F1 are noisy.
- **Simple pathway scoring.** Mean expression is only one option. Rank-based methods (for example AUCell) may capture different information.
- **No gene interactions.** Biological knowledge was used only to build features, not inside the model.
- **Gene selection and cross-validation.** The first grid search did the selection inside the folds (CV macro F1 0.937). For the later models, the genes were selected once on the full training set (0.936). The scores are almost the same, so the effect seems small, and the test set was never used for selection.
- **Random Forest tuning.** The number of trees was chosen using all genes without class weights, while the final model used the reduced gene set with balanced weights.

**Future work:** alternative pathway representations, the contribution of single pathways to the predictions, testing the possible reasons for the performance gap one by one, and models that include pathways or molecular interaction networks directly in their design.

---

## 7. Reproducibility

```bash
pip install scanpy scikit-learn pandas numpy matplotlib seaborn xgboost gseapy
jupyter notebook Single_cell_ML_project_clean.ipynb
```

All random seeds are fixed (`random_state=42`). The dataset is downloaded automatically by `scanpy`. The pathway part needs an internet connection, because `gseapy` downloads the gene sets from Enrichr.

**Tools:** Python, scanpy, scikit-learn, XGBoost, gseapy, pandas, seaborn.
