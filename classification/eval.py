import numpy as np
import scanpy as sc
import pickle
import pandas as pd
from pathlib import Path
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.tree import DecisionTreeClassifier


MYSEED = 2211
np.random.seed(MYSEED)
np.random.RandomState(MYSEED)


sc.settings.figdir = './figures/classifier/'
sc.settings.file_format_figs = 'eps'
sc.settings._vector_friendly = False
sc.settings.autosave = True
sc.settings.autoshow = True
sc.settings._frameon = False


def load_HOM_datasets(dataset_type='variable'):
    E14 = sc.read(Path('ann_data/', 'E14_hom_' + dataset_type + '_genes.h5ad'))
    E13 = sc.read(Path('ann_data/', 'E13_hom_' + dataset_type + '_genes.h5ad'))
    # Convert to dataframes
    E14_df = E14.to_df()
    E13_df = E13.to_df()

    intersection = E14_df.columns.intersection(E13_df.columns)

    # Add Leiden annotations (class) as the class for E14
    if dataset_type == 'all':  # Load annotations from the variable genes dataset
        E14_variable = sc.read(Path('ann_data/', 'E14_hom_variable_genes.h5ad'))
        annotations = E14_variable.obs['leiden_annotations']
    else:
        annotations = E14.obs['leiden_annotations']

    X = E14_df[intersection]
    X_test = E13_df[intersection]

    # Scale X and X_test
    scaler = StandardScaler().fit(X)
    X = scaler.transform(X)
    X_test = scaler.transform(X_test)

    # E14 has more columns than E13. Keep only the intersection of both
    return X, annotations, X_test


def celltypes_percentage(y):
    for celltype in set(y):
        n_pred = len(y)
        n_celltype = len(y[y == celltype])
        print('{}: {}%'.format(celltype, round(n_celltype / n_pred, 2) * 100))


if __name__ == '__main__':
    X, y, X_test = load_HOM_datasets('all')

    # Encode labels
    label_encoder = LabelEncoder().fit(y)
    y_encoded = label_encoder.transform(y)

    # Load CV results from pickle file
    with open('./classification/cv_results.pickle', 'rb') as handle:
        grid_cv = pickle.load(handle)
        cv_results = pd.DataFrame(grid_cv.cv_results_)

    y_pred = grid_cv.predict(X_test)
    y_pred = label_encoder.inverse_transform(y_pred)

    # Percentage of predicted celltypes in E13
    celltypes_percentage(y_pred)

    E13_var = sc.read('ann_data/E13_hom_variable_genes.h5ad')
    E13_var.obs['predictions'] = y_pred
    sc.pl.umap(E13_var, color=['leiden_annotations', 'predictions'], title=['Leiden annotations', 'RandomForest predictions'])

    #
    # Eval with DecisionTrees
    decision_tree_clf = DecisionTreeClassifier(max_depth=6, random_state=MYSEED).fit(X, y)
    y_pred = decision_tree_clf.predict(X_test)
    celltypes_percentage(y_pred)
