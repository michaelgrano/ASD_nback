
import numpy as np
from scipy.io import loadmat, savemat
from sklearn.linear_model import RidgeClassifierCV
import sklearn
from sklearn.model_selection import RepeatedStratifiedKFold
import sklearn.pipeline
import pickle
import argparse
import os
from tqdm import tqdm

def get_classifier(classifier, select_topn_perc=100, n_jobs=1):
    """
    function for getting a classifier with optional feature selection
    args:
        classifier: one of ['ridge','svm','svm_stock'] (svm_stock is non-optimized svm)
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--alt', action='store_true')
    parser.add_argument('--short-irf', action='store_true')
    parser.add_argument('--dtype', type=str, default='hit')
    parser.add_argument('--n-jobs', type=int, default=1)
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--bin-space', type=int, default=1)
    parser.add_argument('--bin-width', type=int, default=10)
    parser.add_argument('--outer-cv', type=int, default=1)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--compute-nulldist', action='store_true')
    args = parser.parse_args()

    # load the data from matlab (see classify.m for how it was saved)
    alt_tag = '_a' if args.alt else ''
    irf_tag = '_2s' if args.short_irf else ''
    data = loadmat(f'data/XY{irf_tag}{alt_tag}_{args.dtype}.mat')
    X = data['data']['X'][0][0]
    Y = data['data']['Y'][0][0].flatten()
    # make equal sized groups given more controls
    X = X[0:2*np.sum(Y==1), :]
    Y = Y[0:2*np.sum(Y==1)]

    ##  all of our parameters
    # sliding window size
    bin_width = args.bin_width
    bin_space = args.bin_space # 1 for all time points
    bin_t0s = np.arange(0, X.shape[1]-bin_width, bin_space) # e.g., np.arange(start, end, space)
    # kfold parameters
    n_repeats = args.outer_cv
    n_splits = 5
    # which classifiers
    classifiers = ['ridge']
    # feature selection for non-sliding window
    topn_perc = 100
    # for null dist
    compute_nulldist = args.compute_nulldist
    null_tag = '_with_nulldist' if compute_nulldist else ''
    n_perms = 10000
    n_perms_use = n_perms if compute_nulldist else 1
    # output data
    fn = f'results/decoding_{args.dtype}{irf_tag}{alt_tag}_{bin_width}-tp-wide_{bin_space}tp-spaced_{n_repeats}outer-cv{null_tag}.pkl'
    fn_matlab = fn.replace('.pkl', '.mat')

    print(f'working on {fn}')

    if os.path.exists(fn) and not args.overwrite:
        with open(fn, 'rb') as f:
            all_data = pickle.load(f)
    else:
        # for splitting the data, let's do 5-fold with 20 repeats
        rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
            random_state=1)

    ###
        full_accs = {}
        for classifier in classifiers:
            full_acc = []
            window_acc = []
            for train_inds, test_inds in rskf.split(X, Y):
                classif, manual_grid_search = get_classifier(classifier, select_topn_perc=topn_perc, n_jobs=args.n_jobs)
                if manual_grid_search:
                    classif = classif.best_estimator_
                classif = classif.fit(X[train_inds,:], Y[train_inds])
                result = classif.predict(X[test_inds,:])
                full_acc.append(np.mean(result == Y[test_inds]))
            full_accs[classifier] = np.mean(full_acc)

        window_accs = {} # dict to store all the classifier sliding accuracies, e.g. classifier_acc['ridge'] is a time-series of accuracies
        for classifier in classifiers:
            window_acc = []
            print(f'beginning classifier {classifier}...')
            for t_0 in tqdm(bin_t0s): #range(X.shape[1]-bin_width):
                acc = []
                for train_inds, test_inds in rskf.split(X, Y):
                    classif, manual_grid_search = get_classifier(classifier, select_topn_perc=topn_perc, n_jobs=args.n_jobs)
                    X_train, Y_train = X[train_inds,t_0:t_0+bin_width], Y[train_inds]
                    Y_train = np.random.permutation(Y_train)
                    classif = classif.fit(X[train_inds,t_0:t_0+bin_width], Y[train_inds])
                    if manual_grid_search:
                        classif = classif.best_estimator_
                    result = classif.predict(X[test_inds,t_0:t_0+bin_width])
                    acc.append(np.mean(result == Y[test_inds]))
                window_acc.append(np.mean(acc))
            window_accs[classifier] = np.array(window_acc)

        if compute_nulldist:
            full_accs_null = {}
            for classifier in classifiers:
                full_accs_null[classifier] = []
                print(f'beginning classifier {classifier} null permutations...')
                for nullperm in tqdm(range(n_perms)):
                    full_acc = []
                    window_acc = []
                    for train_inds, test_inds in rskf.split(X, Y):
                        classif, manual_grid_search = get_classifier(classifier, select_topn_perc=topn_perc, n_jobs=args.n_jobs)
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

            window_accs_null = {} # dict to store all the classifier sliding accuracies, e.g. classifier_acc['ridge'] is a time-series of accuracies
            for classifier in classifiers:
                window_accs_null[classifier] = []
                print(f'beginning classifier {classifier} temporal null permutations...')
                for nullperm in tqdm(range(n_perms)):
                    window_acc = []
                    for t_0 in bin_t0s:
                        acc = []
                        for train_inds, test_inds in rskf.split(X, Y):
                            classif, manual_grid_search = get_classifier(classifier, select_topn_perc=topn_perc, n_jobs=args.n_jobs)
                            # if we were doing null dist., we would permute data
                            X_train, Y_train = X[train_inds,t_0:t_0+bin_width], Y[train_inds]
                            X_test, Y_test = X[test_inds,t_0:t_0+bin_width], Y[test_inds]
                            if compute_nulldist:
                                Y_train = np.random.permutation(Y_train)
                            classif = classif.fit(X_train, Y_train)
                            if manual_grid_search:
                                classif = classif.best_estimator_
                            result = classif.predict(X_test)
                            acc.append(np.mean(result == Y_test))
                        window_acc.append(np.mean(acc))
                    window_accs_null[classifier].append(np.array(window_acc))
                window_accs_null[classifier] = np.array(window_accs_null[classifier])

                all_data = dict(full_accs=full_accs, window_accs=window_accs,
                        full_accs_null=full_accs_null, window_accs_null=window_accs_null)
        else:
            all_data = dict(full_accs=full_accs, window_accs=window_accs)
        with open(fn, 'wb') as f:
            pickle.dump(all_data, f)

        savemat(fn_matlab, all_data)

    if args.plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.set_xlabel('Beginning of 10 time point bin')
        for classifier in classifiers:
            ax.plot(all_data['window_accs'][classifier], label=classifier)
        ax.legend()
        ax.set_ylabel('Accuracy')
        fig.savefig(f'figures/classify_autism{alt_tag}_{args.dtype}.png', bbox_inches='tight')

# all_data is a grand dictionary containing full_accs and window_accs
# full_accs stores the accuracies computed across the full time series for each classifier
# window_accs stores the sliding window accuracies for each classifier
