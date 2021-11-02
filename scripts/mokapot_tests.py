import mokapot
import pandas as pd
import os
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

np.random.seed(42)

pin_file = "results_SA_RTdiff_absRT.pin"
database = "iPRG2015.TargDecoy.fasta"


# Analyze PSMs from a single file with mokapot using default SVM model
"""
# Read the PSMs from the PepXML file:
psms = mokapot.read_pin(pin_file)

#psms.add_proteins(database)

# Conduct the mokapot analysis:
results, models = mokapot.brew(psms)

# Save the results to two tab-delimited files
# "mokapot.psms.txt" and "mokapot.peptides.txt"
result_files = results.to_txt()
"""

# Random Forest model used by AlphaPept (Strauss et al.)
"""
scaler = StandardScaler
rfc = RandomForestClassifier()  # class_weight={False:1,True:5},
# Initiate scaling + classification pipeline
pipeline = Pipeline([('scaler', scaler), ('clf', rfc)])
parameters = {'clf__max_depth': [5,25,50], 'clf__max_leaf_nodes': [150,200,250]}
# Setup grid search framework for parameter selection and internal cross validation
cv = GridSearchCV(pipeline, param_grid=parameters, cv=5, scoring="accuracy",
                  verbose=0, return_train_score=True, n_jobs=-1)
"""


# Class weight param grid for the default SVM
grid = {"class_weight": [{0: neg, 1: pos} for neg in (0.1, 1, 10) for pos in (0.1, 1, 10)]}
#grid = {"max_depth": [5,25,50], "max_leaf_nodes": [150,200,250]}

# Model class of reimplemented Percolator model. Modify to test different models provided by sklearn package
class TestModel(mokapot.model.Model):
    def __init__(
        self,
        scaler=None,
        train_fdr=0.01,
        max_iter=10,
        direction=None,
        override=False,
        subset_max_train=None,
    ):
        model = RandomForestClassifier()
        estimator = GridSearchCV(model, param_grid=grid, refit=False, cv=5)
        #estimator = cv

        # Default model used by mokapot
        #svm_model = LinearSVC(dual=False)
        #estimator = GridSearchCV(
        #    svm_model, param_grid=GRID, refit=False, cv=3
        #)

        super().__init__(
            estimator=estimator,
            scaler=scaler,
            train_fdr=train_fdr,
            max_iter=max_iter,
            direction=direction,
            override=override,
            subset_max_train=subset_max_train,
        )

# Analyze PSMs from a single file with protein-level results

# Read the PSMs from a pin file:
psms = mokapot.read_pin(pin_file)

# Provide the protein sequences, set decoy prefix:
psms.add_proteins(database, decoy_prefix="DECOY_sp")

# Conduct the mokapot analysis:
results, models = mokapot.brew(psms, model=TestModel())

# Save the results to three tab-delimited files
# "mokapot.psms.txt", "mokapot.peptides.txt", and "mokapot.proteins.txt"
result_files = results.to_txt()

# Print out number of hits showing q-value <= 0.01 (1% FDR)
moka_psms = pd.read_table("mokapot.psms.txt")

moka_passing = (moka_psms["mokapot q-value"] <= 0.01).sum()
print(moka_passing)