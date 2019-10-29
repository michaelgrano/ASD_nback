
import numpy as np
from scipy.io import loadmat, savemat
from sklearn.linear_model import RidgeClassifierCV
import sklearn
from sklearn.model_selection import RepeatedStratifiedKFold
import scipy
import sklearn.pipeline
import pickle
import argparse
import os
from tqdm import tqdm

def get_classifier(classifier, select_topn_perc=100, n_jobs=1):
    """
    function for getting a classifier with optional feature selection
    args:
        classifier: one of ['lr' 'ridge','svm','svm_stock'] (svm_stock is non-optimized svm)
        select_topn_perc: top n percentile of features to select with F-test
        n_jobs: number of parallel jobs to use
    returns:
        clf: the classifier
        manual_grid_search: whether we need to do manual grid search
    """
    manual_grid_search = False
    selector = None
    if classifier == 'ridge':
        clf_args = {'alphas': np.logspace(-3, 6, num=100), 'cv': None}
        estimator = sklearn.linear_model.RidgeClassifierCV(**clf_args)
        selector = sklearn.feature_selection.SelectPercentile(percentile=select_topn_perc)
        clf = sklearn.pipeline.Pipeline([('selector', selector),
                    ('ridge', estimator)])
    elif classifier == 'svm':
        estimator = sklearn.svm.SVC(kernel='linear')
        selector = sklearn.feature_selection.SelectPercentile(percentile=select_topn_perc)
        pipeline = sklearn.pipeline.Pipeline([('selector', selector),
                    ('svc', estimator)])
        param_grid = {'svc__C': np.logspace(-3, 6, num=10)}
        clf = sklearn.model_selection.GridSearchCV(pipeline, param_grid, cv=2, n_jobs=n_jobs)
        manual_grid_search = True
    elif classifier == 'svm_stock':
        estimator = sklearn.svm.SVC(kernel='linear')
        selector = sklearn.feature_selection.SelectPercentile(percentile=select_topn_perc)
        clf = sklearn.pipeline.Pipeline([('selector', selector),
                    ('svc', estimator)])
    elif classifier == 'lr':
        clf_args = {'penalty':'none', 'solver':'saga', 'class_weight': 'balanced'}
        estimator = sklearn.linear_model.LogisticRegression(**clf_args)
        selector = sklearn.feature_selection.SelectPercentile(percentile=select_topn_perc)
        clf = sklearn.pipeline.Pipeline([('selector', selector),
                    ('lr', estimator)])
    return clf, manual_grid_search

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--distractors', action='store_true', help= 'use to run distractors condition')
    parser.add_argument('--distractor-diff', action='store_true', help='use to take metric difference b/w distractors absent/present conditions')
    parser.add_argument('--dtype', type=str, default='hit', choices=['hit', 'FA', 'miss'], help='trial type')
    parser.add_argument('--n-jobs', type=int, default=1, help='number of parallel workers to use per classifier. usually not much speedup for >1')
    parser.add_argument('--outer-cv', type=int, default=20, help='# of outer loops in cross-validation')
    parser.add_argument('--k-folds', type=int, default=5, help='# of folds in inner cross-validation loop')
    parser.add_argument('--measure', type=str, choices=['MAD', 'MaxAD', 'MeanAD', 'MADs'], default='MAD', help='IRF metric for decoding on; MADs uses all three')
    parser.add_argument('--overwrite', action='store_true', help=' ')
    parser.add_argument('--no-nulldist', action='store_true', help='use to disable nulldist compuutation and save a lot of time')
    args = parser.parse_args()

    # HARDCODED params
    # kfold parameters
    n_repeats = args.outer_cv
    k_folds = args.k_folds
    classifiers = ['lr'] # which classifiers
    # for null dist
    compute_nulldist = not args.no_nulldist
    null_tag = '_with_nulldist' if compute_nulldist else ''
    n_perms = 10000

    # load the data from matlab (see classify.m for how it was saved)
    alt_tag = '_a' if args.distractors else '_distractors-diff' if args.distractor_diff else ''
    data = loadmat(f'data/XY{alt_tag}_{args.dtype}.mat')
    X = data['data']['X'][0][0]
    Y = data['data']['Y'][0][0].flatten()

    # output data
    fn = f'results/decoding_{args.dtype}{alt_tag}_{args.measure}_{n_repeats}outer-cv{null_tag}.pkl'
    fn_matlab = fn.replace('.pkl', '.mat')
    print(f'working on {fn}...')

    # compute MAD
    if args.distractor_diff:
        # we calculated MAD ahead of time, just make 2D
        if args.measure != 'MAD':
            raise ValueError()
    else:
        if args.measure == 'MAD':
            X = np.expand_dims(np.median(np.abs(X - np.expand_dims(np.median(X, axis=1),axis=1)), axis=1),axis=1)
        elif args.measure == 'MaxAD':
            X = np.expand_dims(np.abs(np.max(X, axis=1) - np.min(X, axis=1)),1)
        elif args.measure == 'MeanAD':
            X = np.expand_dims(np.mean(np.abs(X - np.expand_dims(np.mean(X, axis=1),axis=1)), axis=1),axis=1)
        elif args.measure == 'MADs':
            X = np.stack((np.median(np.abs(X - np.expand_dims(np.median(X, axis=1),axis=1)), axis=1),
                                np.abs(np.max(X, axis=1) - np.min(X, axis=1)),
                                np.mean(np.abs(X - np.expand_dims(np.mean(X, axis=1),axis=1)), axis=1)),
                                axis=1)

    if os.path.exists(fn) and not args.overwrite:
        with open(fn, 'rb') as f:
            all_data = pickle.load(f)
    else:
        # for splitting the data, let's do 5-fold with 20 repeats
        rskf = RepeatedStratifiedKFold(n_splits=k_folds, n_repeats=n_repeats, random_state=1)

        full_accs = {}
        for classifier in classifiers:
            full_acc = []
            window_acc = []
            for train_inds, test_inds in rskf.split(X, Y):
                classif, manual_grid_search = get_classifier(classifier, n_jobs=args.n_jobs)
                if manual_grid_search:
                    classif = classif.best_estimator_
                classif = classif.fit(X[train_inds,:], Y[train_inds])
                result = classif.predict(X[test_inds,:])
                full_acc.append(np.mean(result == Y[test_inds]))
            full_accs[classifier] = np.mean(full_acc)

        if compute_nulldist:
            full_accs_null = {}
            for classifier in classifiers:
                full_accs_null[classifier] = []
                print(f'beginning classifier {classifier} null permutations...')
                for nullperm in tqdm(range(n_perms)):
                    full_acc = []
                    for train_inds, test_inds in rskf.split(X, Y):
                        classif, manual_grid_search = get_classifier(classifier, n_jobs=args.n_jobs)
                        if manual_grid_search:
                            classif = classif.best_estimator_
                        X_train, Y_train = X[train_inds,:], Y[train_inds]
                        X_test, Y_test = X[test_inds,:], Y[test_inds]
                        if compute_nulldist:
                            Y_train = np.random.permutation(Y_train)
                        classif = classif.fit(X_train, Y_train)
                        result = classif.predict(X_test)
                        full_acc.append(np.mean(result == Y_test))
                    full_accs_null[classifier].append(np.mean(full_acc))
                full_accs_null[classifier] = np.array(full_accs_null[classifier])

            all_data = dict(full_accs=full_accs, full_accs_null=full_accs_null)
        else:
            all_data = dict(full_accs=full_accs)
        with open(fn, 'wb') as f:
            pickle.dump(all_data, f)

        savemat(fn_matlab, all_data)
