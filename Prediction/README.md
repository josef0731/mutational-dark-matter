# Predictors of variant effects trained on DMS data

This folder contains data and code on predictors of variant effects trained over (experimental/computational) DMS data. **Please read [here](https://github.com/josef0731/mutational-dark-matter/blob/main/Prediction/using-predictors.ipynb) for a Quick-Start introduction on how to apply the predictors to data.**

## Background

Briefly, we have variant effect predictors trained on two datasets:,

1. `DMSexp`: trained on experimental DMS data covering 10 proteins.
2. `EVmutation`: trained on computational DMS data given by a scoring system implemented in [EVmutation](https://marks.hms.harvard.edu/evmutation/).

For each, we have trained predictors using these different combinations of features:

* AA + Cons (i.e. Conservation) + MutSig + QSASA (i.e. Solvent accessibility)
* AA + Cons + MutSig
* AA + MutSig + QSASA
* AA + MutSig
* MutSig only
* AA only

Additionally, for `DMSexp` and `EVmutation`, we also trained a version of predictors where variants residing in the protein core were removed from training, to evaluate the importance of protein core positions to learning variant effects. We call these the *no-core* predictors.

All predictors were trained using the Gradient Boosting Classifier method in scikit-learn, and can be loaded from the pickle objects included in the repository - `using-predictors.ipynb` is a Jupyter notebook with examples on how to annotate variant data to apply the predictors.

## Description of the files

### Folder DMSexp

Predictors trained on experimental DMS data.

`download_MAVEdb_final.ipynb`: Jupyter notebook with code to interrogate MAVEdb and download data from the database.

`DMSexp_prepare_data_final.ipynb`: Jupyter notebook with code to prepare, clean and annotate the variants for training/testing the predictors.

`DMSexp_prediction-GBM_final.ipynb`: Jupyter notebook containing details on training and testing predictors trained on the experimental DMS data.

Subfolder `setup`: with separate subfolders containing python `pickle` object of the Gradient Boosting Classifier (scikit-learn objects) to be loaded and applied on variant datasets. Each setup (i.e. tested combination of features) contain:

1. `GBM.pickle` - i.e. the model trained on all available training set variants.
2. `GBM-no-core.pickle` - i.e. the model trained on the training set with variants residing in the protein core removed.

### Folder EVmutation

Predictors trained on computational DMS data generated using the EVmutation scoring method. Directly downloaded from EVmutation website.

`EVmutation_prepare_data_final.ipynb`: Jupyter notebook with code to prepare, clean and annotate EVmutation variants for training/testing the predictors.

`EVmutation_prediction-GBM_final.ipynb`: Jupyter notebook on training and testing predictors using the EVmutation scores.

Subfoler `setup`: with separate subfolders containing python `pickle` object of the Gradient Boosting Classifier (scikit-learn objects) to be loaded and applied on variant datasets. Each setup (i.e. tested combination of features) contain:

1. `GBM.pickle` - i.e. the model trained on all available training set variants.
2. `GBM-no-core.pickle` - i.e. the model trained on the training set with variants residing in the protein core removed.
3. `GBM-excludeBRCA1-PTEN.pickle` - i.e. the model trained on the training set with variants in the BRCA1 and PTEN proteins removed from training (used for testing the predictors against experimental DMS data on BRCA1 and PTEN).
4. `GBM-excludeBRCA1-PTEN.pickle` - same as (3) but further removing those variants not in the two named proteins, but are residing in the protein core.

### Test case

Here as illustration we provide examples on applying the predictors against the protein BRCA1. The results are also discussed in the manuscript.

`Findlay_et_al_BRCA_DMS.txt`: Experimental DMS measurements obtained from [Findlay et al Nature].

`BRCA1_EVmutation.tsv`: BRCA1 variant scores from EVmutation.

`BRCA1_ZoomVar.tsv`: BRCA1 protein structural mapping obtained from ZoomVar.

`BRCA1_DMS_annotated.csv`: Annotated BRCA1 variants ready for applying the predictors.

`BRCA1_DMS_annotated.pickle`: Same data as the CSV but in pickle format. Converted into a (X, Y) duple where X is the feature matrix (in numpy array format) and Y is the variant effect binary label.

`using-predictors.ipynb`: **Jupyter notebook containing examples of applying the predictors. Read this to learn how to use our predictors!**
 
`helper.py`: helper functions used in the `using-predictors.ipynb` notebook.
