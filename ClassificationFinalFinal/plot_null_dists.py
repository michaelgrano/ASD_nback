
import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse

def plot_null_dists(bin_space, bin_width, cv_repeats, short_irf):
    if short_irf:
        xtime = np.arange(1/40, 2-bin_width/40 +10e-4, bin_space/40)
    else:
        xtime = np.arange(1/40, 4-bin_width/40 +10e-4, bin_space/40)

    bonf = int(np.floor(10000*(.05/len(xtime))))

    cv_tag = '' if cv_repeats is 1 else f'_{cv_repeats}outer-cv'
    irf_tag = '_2s' if short_irf else ''
    fig, axs = plt.subplots(2,4)

    print('\n \n ----------------------------------------------')
    if short_irf:
        print('2s IRF')
    else:
        print('4s IRF')
    for ii, dtype in enumerate(['hit', 'miss', 'FA']):

        with open(f'results/decoding_{dtype}{irf_tag}_a_{bin_width}-tp-wide_{bin_space}tp-spaced{cv_tag}_with_nulldist.pkl', 'rb') as f:
            res_a = pickle.load(f)

        # just print to console non-sliding window results
        if ii==0:
            nulldist = np.sort(res_a['full_accs_null']['ridge'], axis=0)
            pval = (10000-np.argmin(abs(nulldist - res_a['full_accs']['ridge'])))/10000
            print(f"{dtype} \n \t distractors \n \t acc={res_a['full_accs']['ridge']:.03f}, p={pval:.03f}")

        # ax = axs[1,ii]
        fig, ax = plt.subplots(figsize=(4.2,3.2))
        nulldist_a = np.sort(res_a['window_accs_null']['ridge'], axis=0)
        ax.plot(xtime, res_a['window_accs']['ridge'], label='ridge classifier')
        ax.plot(xtime, nulldist_a[10000-bonf,:], 'r--', label='p<0.05 bonferroni')
        ax.plot(xtime, nulldist_a[9500,:], 'r-.', label='p<0.05 uncorrected')
        ax.plot(xtime, nulldist_a[5000,:], 'r-', label='mean null dist')
        ax.set_ylim([.35,.75])
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Classification accuracy')
        ax.set_title(f'Distractors: {dtype}')
        ax.legend(fontsize='x-small')
        fig.savefig(f'figures/decoding_{dtype}{irf_tag}_a_{bin_width}-tp-wide_{bin_space}tp-spaced{cv_tag}_with_nulldist.png', dpi=200, bbox_inches='tight')
        # plt.show()

        with open(f'results/decoding_{dtype}{irf_tag}_{bin_width}-tp-wide_{bin_space}tp-spaced{cv_tag}_with_nulldist.pkl', 'rb') as f:
            res = pickle.load(f)

        # just print to console non-sliding window results
        nulldist = np.sort(res['full_accs_null']['ridge'], axis=0)
        pval = (10000-np.argmin(abs(nulldist - res['full_accs']['ridge'])))/10000
        print(f"\t no distractors \n \t acc={res['full_accs']['ridge']:.03f}, p={pval:.03f}")

        # ax = axs[0,ii]
        fig, ax = plt.subplots(figsize=(4.2,3.2))
        nulldist = np.sort(res['window_accs_null']['ridge'], axis=0)
        ax.plot(xtime, res['window_accs']['ridge'], label='ridge classifier')
        ax.plot(xtime, nulldist[10000-bonf,:], 'r--', label='p<0.05 bonferroni')
        ax.plot(xtime, nulldist[9500,:], 'r-.', label='p<0.05 uncorrected')
        ax.plot(xtime, nulldist[5000,:], 'r-', label='mean null dist')
        ax.set_ylim([.35,.75])
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Classification accuracy')
        ax.set_title(f'No distractors: {dtype}')
        # ax.legend(fontsize='x-small')
        fig.savefig(f'figures/decoding_{dtype}{irf_tag}_{bin_width}-tp-wide_{bin_space}tp-spaced{cv_tag}_with_nulldist.png', dpi=200, bbox_inches='tight')
        # plt.show()

    # fig.legend(fontsize='small')
    # fig.savefig(f'figures/decoding_all_{bin_width}-tp-wide_{bin_space}tp-spaced{cv_tag}.png', dpi=200, bbox_inches='tight')
    # plt.show()


    for ii, dtype in enumerate(['hit', 'miss', 'FA']):

        with open(f'results/decoding_{dtype}{irf_tag}_{bin_width}-tp-wide_{bin_space}tp-spaced{cv_tag}_with_nulldist.pkl', 'rb') as f:
            res = pickle.load(f)
        with open(f'results/decoding_{dtype}{irf_tag}_a_{bin_width}-tp-wide_{bin_space}tp-spaced{cv_tag}_with_nulldist.pkl', 'rb') as f:
            res_a = pickle.load(f)

        # ax = axs[0,ii]
        fig, ax = plt.subplots()
        nulldist = np.sort(res_a['window_accs_null']['ridge'] - res['window_accs_null']['ridge'], axis=0)
        ax.plot(xtime, res_a['window_accs']['ridge'] - res['window_accs']['ridge'], label='ridge classifier')
        ax.plot(xtime, nulldist[10000-bonf,:], 'r--', label='p<0.05 bonferroni')
        ax.plot(xtime, nulldist[9750,:], 'r-.', label='p<0.05 uncorrected')
        ax.plot(xtime, nulldist[5000,:], 'r-', label='mean null dist')
        # ax.legend(fontsize='x-small')
        ax.plot(xtime, nulldist[250,:], 'r--', label='p<0.05 uncorrected')
        ax.plot(xtime, nulldist[bonf,:], 'r-.', label='p<0.05 bonferroni')
        # ax.set_ylim([.35,.75])
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Classification accuracy')
        ax.set_title(f'(Distractor) - (no distractors): {dtype}')
        fig.savefig(f'figures/group_diff_decoding_{dtype}{irf_tag}_{bin_width}-tp-wide_{bin_space}tp-spaced{cv_tag}_with_nulldist.png', dpi=200, bbox_inches='tight')
        # plt.show()
    print('\n \n ----------------------------------------------')


def print_MAD_results(MAD_measure, cv_repeats=20, short_irf=False):
    print(f'--------------{MAD_measure}---------------')
    irf_tag = '_2s' if short_irf else ''
    res = {}
    print('\n------------------')
    for dtype in ['hit', 'miss', 'FA']:
        print(f'{dtype} \n')
        for distractors in ['', '_a']:
            dtag = 'distractors' if 'a' in distractors else 'no distractors'
            with open(f'results/decoding_{dtype}{irf_tag}{distractors}_{MAD_measure}_{cv_repeats}outer-cv_with_nulldist.pkl', 'rb') as f:
                #[f'{dtype}{distractors}']
                res = pickle.load(f)
                acc = res['full_accs']['lr']
                nulldist = np.sort(res['full_accs_null']['lr'], axis=0)
                pval = (10000-np.argmin(abs(nulldist - acc)))/10000
                print(f"{dtag}:\n acc={acc:.04f},p={pval:.04f}")
        print('------------------')


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

def plot_full_window_results(bin_space, bin_width, cv_repeats, short_irf):
    irf_tag = '_2s' if short_irf else ''
    cv_tag = '' if cv_repeats is 1 else f'_{cv_repeats}outer-cv'
    row_labels = ['hit', 'miss', 'FA']
    col_labels= ['no distractors', 'distractors']
    datatxt = []
    for dtype in ['hit', 'miss', 'FA']:
        rowtxt = []
        for distractors in ['', '_a']:
            dtag = 'distractors' if 'a' in distractors else 'no distractors'
            with open(f'results/decoding_{dtype}{irf_tag}{distractors}_{bin_width}-tp-wide_{bin_space}tp-spaced{cv_tag}_with_nulldist.pkl', 'rb') as f:
                res = pickle.load(f)
                acc = res['full_accs']['ridge']
                nulldist = np.sort(res['full_accs_null']['ridge'], axis=0)
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
    fig.savefig('figures/full_window_table.png', bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bin-space', type=int, default=10)
    parser.add_argument('--bin-width', type=int, default=20)
    parser.add_argument('--cv_repeats', type=int, default=3)
    parser.add_argument('--short-irf', action='store_true')
    parser.add_argument('--MAD-measure', type=str, default='MAD', choices=['MAD', 'MaxAD', 'MeanAD', 'MADs'])
    parser.add_argument('--plot-all', action='store_true')
    args = parser.parse_args()

    if args.plot_all:
        bin_space = 10
        bin_widths = [5,10,20,40]
        cv_repeats = 20
        short_irf = False
        # for bin_width in bin_widths:
        #     plot_null_dists(bin_space, bin_width, cv_repeats, short_irf)
        plot_MAD_results(cv_repeats, short_irf)
        # plot_full_window_results(10, 20, cv_repeats, short_irf)
    else:
        plot_null_dists(args.bin_space, args.bin_width, args.cv_repeats, args.short_irf)
        print_MAD_results(args.MAD_measure, args.cv_repeats, args.short_irf)
