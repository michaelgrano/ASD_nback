
import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def plot_MAD_results(cv_repeats=20, short_irf=False, MAD_measure='MAD'):
    irf_tag = '_2s' if short_irf else ''
    row_labels = ['hit', 'hit-distractors', 'hit-distractor-diff', 'FA', 'FA-distractors', 'FA-distractor-diff', 'miss', 'miss-distractors', 'miss-distractor-diff']
    row_labels = ['hit', 'FA', 'miss']
    col_labels = [r'MAD$_{distractors-absent}$', r'MAD$_{distractors-present}$', r'$\Delta_{distractors}$MAD']
    # col_labels= ['MAD'] #['MedAD', 'MeanAD', 'MaxAD', 'MADs']
    datatxt = []
    for dtype in ['hit', 'FA', 'miss']:
        rowtxt = []
        for distractors in ['', '_a', '_distractors-diff']:
            with open(f'results/decoding_{dtype}{irf_tag}{distractors}_{MAD_measure}_{cv_repeats}outer-cv_with_nulldist.pkl', 'rb') as f:
                #[f'{dtype}{distractors}']
                res = pickle.load(f)
                acc = res['full_accs']['lr']
                nulldist = np.sort(res['full_accs_null']['lr'], axis=0)
                pval = (10000-np.argmin(abs(nulldist - acc)))/10000
                if acc > np.max(nulldist):
                    rowtxt.append(r'$\bf\mu={:.03f}$, $p<0.0001$'.format(acc, pval))
                elif pval < (.05/24):
                    rowtxt.append(r'$\bf\mu={:.03f}$, $p={:.04f}$'.format(acc, pval))
                else:
                    rowtxt.append(r'$\mu={:.03f}$, $p={:.04f}$'.format(acc, pval))

        datatxt.append(rowtxt)
    fig, ax = plt.subplots(figsize=(6,3))
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.axis('off')
    ax.table(cellText=datatxt, rowLabels=row_labels, colLabels=col_labels, loc='center')
    plt.show()
    fig.savefig(f'figures/{MAD_measure}_table.png', bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cv_repeats', type=int, default=20)
    parser.add_argument('--MAD-measure', type=str, default='MAD', choices=['MAD', 'MaxAD', 'MeanAD', 'MADs'])
    args = parser.parse_args()
    opt = args.__dict__

    os.makedirs('figures', exist_ok=True)
    plot_MAD_results(**opt)
