import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pickle
from pathlib import Path
from sklearn import svm
from sklearn.model_selection import validation_curve, GridSearchCV
from sklearn.metrics import classification_report, plot_confusion_matrix
from sklearn.pipeline import Pipeline
from sklearn.linear_model import SGDClassifier, LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.tree import DecisionTreeClassifier

MYSEED = 2211
np.random.seed(MYSEED)
np.random.RandomState(MYSEED)


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


def print_classif_report(y_true, y_pred, labels):
    print(classification_report(y_true, y_pred))


def confusion_matrix(classifier, X_test, y_test, class_names, classifier_name):
    np.set_printoptions(precision=2)

    disp = plot_confusion_matrix(classifier, X_test, y_test, display_labels=class_names, cmap=plt.cm.Blues,
                                 normalize='true', xticks_rotation=20)
    disp.ax_.set_title('Normalized confusion matrix - ' + classifier_name)
    print(disp.confusion_matrix)

    plt.tight_layout()
    plt.show()


def train_val_score(param_range, X_train, y_train):
    train_scores, test_scores = validation_curve(
        svm.SVC(), X_train, y_train, param_name='gamma', param_range=param_range,
        scoring='accuracy', n_jobs=1)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)

    plt.title('Validation Curve with SVM')
    plt.xlabel(r'$\gamma$')
    plt.ylabel('Score')
    plt.ylim(0.0, 1.1)
    lw = 2
    plt.semilogx(param_range, train_scores_mean, label='Training score',
                 color='darkorange', lw=lw)
    plt.fill_between(param_range, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.2,
                     color='darkorange', lw=lw)
    plt.semilogx(param_range, test_scores_mean, label='Cross-validation score',
                 color='navy', lw=lw)
    plt.fill_between(param_range, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.2,
                     color='navy', lw=lw)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()


def celltypes_percentage(y):
    for celltype in set(y):
        n_pred = len(y)
        n_celltype = len(y[y == celltype])
        print('{}: {}%'.format(celltype, round(n_celltype / n_pred, 2) * 100))


if __name__ == '__main__':
    # X, X_test, annotations = load_HOM_datasets('all')
    # X_train, X_val, y_train, y_val = train_test_split(X, annotations, test_size=0.2, random_state=MYSEED)
    X, y, _ = load_HOM_datasets('all')
    celltypes_percentage(y)

    # Encode labels
    label_encoder = LabelEncoder().fit(y)
    y_encoded = label_encoder.transform(y)

    # ==============================
    #       CLASSIFIERS
    # ==============================

    pipeline = Pipeline([
        # ('scal', StandardScaler()),
        ('clf', svm.SVC()),
    ], verbose=True)

    # Grid search
    param_grid = [
        {
            'clf': (SGDClassifier(penalty='l2', alpha=1, random_state=MYSEED),),
        }, {
            'clf': (RandomForestClassifier(n_estimators=200, random_state=MYSEED),),
        }, {
            'clf': (DecisionTreeClassifier(max_depth=6, random_state=MYSEED),),
            'clf__criterion': ('gini', 'entropy'),
            # 'clf__max_depth': (None, 4, 6, 8, 12),
            # 'clf__min_samples_leaf': (1, 5, 10, 15, 20),
            'clf__class_weight': (None, 'balanced')
        }, {
            'clf': (LogisticRegression(max_iter=200, C=1, random_state=MYSEED),),
        },
    ]
    grid = GridSearchCV(pipeline, param_grid=param_grid, n_jobs=20, scoring=['f1_macro', 'accuracy'],
                        refit='f1_macro', verbose=3, cv=3)
    grid.fit(X, y_encoded)
    print('The best model is ', grid.best_estimator_)
    print("The best parameters are %s with a score of %0.2f" % (grid.best_params_, grid.best_score_))
    # The best parameters are {'C': 0.1, 'gamma': 0.1} with a score of 0.35
    # y_pred_1 = grid.predict(X)
    #

    # svm_rbf = svm.SVC(C=0.1, gamma=0.1).fit(X_train, y_train)
    # y_pred_1 = svm_rbf.predict(X_val)
    # confusion_matrix(svm_rbf, X_val, y_val, class_labels, classifier_name='SVM (rbf)')
    # train_val_score([0.001, 0.01, 0.1], X_train, y_train)

    # print_classif_report(y_val, y_pred_1, class_labels)

    # Write CV results to pickle file
    with open('./classification/cv_results.pickle', 'wb') as handle:
        pickle.dump(grid, handle, protocol=pickle.HIGHEST_PROTOCOL)
