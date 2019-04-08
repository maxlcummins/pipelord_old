from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram

import numpy as np
import pandas as pd
import glob

import matplotlib
matplotlib.use('Agg')   #avoid problems with headless environments

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import gridspec


def binning(data, step):
    """
    Bin an array by step width, report the median of each bin
    :param data: array of ordred values to bin
    :param step: the bin step size
    :return: medians, bins
    """
    meds = []
    c = []
    n = 0
    bins = []
    for ir in xrange(data.shape[0]):
        if data[ir, 0] > (n+1)*step:
            if len(c) > 0:
                meds.append(np.median(c))
            else:
                meds.append(0)
            bins.append((n+1)*step)
            n += 1
            c = []
        else:
            c.append(data[ir, 1])
    return np.array(meds), np.array(bins)


def get_block_medians(glob_path, step_size):
    """
    Read in the per-position coverage data for each sample and bin to a corser size

    :param glob_path: path which to scan for files
    :param step_size: the binning size in bp
    :return: pandas dataframe of binned depths (median values), bins for each set of medians
    """

    meds = []
    bins = []
    smpl_names = []

    # manually parse each organism directory. Manually meaning, change this path
    # for each one. Could be made a loop, but unnecessary for 3 runs
    for n, fn in enumerate(glob.glob('{}/*_coverage.txt'.format(glob_path))):

        print '\t {} found {}...'.format(n+1, fn)

        # read the data file
        smpl = fn.split('/')[-1][:-len('_coverage.txt')]
        smpl_names.append(smpl)

        data = pd.read_csv(fn, delimiter='\t', names=['_', 'x', 'depth'])
        data = data.as_matrix()[:, 1:]
        for row in data:
            if row[1] >= 10 and row[1] < 20:
                row[1] = 0.5
            elif row[1] >= 20:
                row[1] = 1
        m, b = binning(data, step_size)
        meds.append(m)
        bins.append(b)

    max_bins = 0
    for bi in bins:
        if max(bi) > max_bins:
            max_bins = max(bi)
    
    for i in xrange(len(bins)):
        if max(bins[i]) < max_bins:
            print 'Warning, data file {} had fewer bins that the largest. Extra empty bins added'.format(i+1)
            bi = bins[i].tolist()
            mi = meds[i].tolist()
            while max(bi) < max_bins:
                bi.append(bi[-1] + step_size)
                mi.append(0)
            bins[i] = np.array(bi)
            meds[i] = np.array(mi)
        

    return pd.DataFrame(np.vstack(meds).T, columns=smpl_names), bins


def run_hclust(outname, meds, bins, step_size, tick_spc, use_polo=True, save_plot=False):
    """
    Cluster and plot the binned data matrix. This calls optimal leaf ordering
    algorithm (polo) by default, which has significant time-complexity.

    :param outname: pdf output file
    :param meds: median values
    :param bins: bins for medians
    :param step_size: size of a step
    :param tick_spc: tick spacing in "# bins"
    :param olo: reorder tips of tree with polo
    :param savePlot: save plot to file, otherwise plot to screen
    """

    # just for indexing
    names = meds.columns

    # normalise rows by their median value
    D = meds.as_matrix().T

#    D = (D.T/np.median(D, 1)).T
    # center rows on their medians
#    D = (D.T - np.median(D, 1)).T
    

    # clustering options. We will use correlation as measure of distance
    # and complete linkage clustering
    metric = 'euclidean'
    method = 'complete'

    # calculate
    Y = linkage(D, method=method, metric=metric)

    # additionally, find optimal leaf ordering
    if use_polo:
        import polo
        print '\tcalculating optimal leaf ordering...'
        Y = polo.optimal_leaf_ordering(Y, pdist(D, metric=metric))

    # now we do some plotting
    fig = plt.figure(figsize=(12, 0.25*len(meds.columns)), dpi=150)
    gs = gridspec.GridSpec(2, 2, width_ratios=[7, 1], height_ratios= [0.2,10])

    axmatrix = plt.subplot(gs[2])
    axmatrix.set_xlabel('genomic coord (kbp)')
    axcolor = plt.subplot(gs[0])
    axcolor.set_title('log median relative abundance ')
    axdend = plt.subplot(gs[3])

    # calculate and plot the dendrogram
    Z = dendrogram(Y, ax=axdend, orientation='right', no_labels=True, color_threshold=0 )

    # the tips (leaves) of the tree become the order for rows
    idx = Z['leaves']

    # reorder rows
    D = np.log(D[idx, :]+1)

    # the largest value in the matrix will set the upper and lower bounds
    # for the heatmap color-range. This assures 0 as the center.
    vmin = np.percentile(D, 2)
    vmax = np.percentile(D, 98)

    # plot the matrix, 5% extra bounds above and below on colour range
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap='PuRd',
                          norm=colors.Normalize(vmin=vmin, vmax=vmax))

    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # try and get some useful axis labels
    #xticks = np.linspace(0, max(bins) - max(bins) % step_size, 5)
    #print xticks

    print '\tticks will be every {}x{} = {} bp'.format(step_size, tick_spc, tick_spc*step_size)
    xticks = np.arange(0, len(bins), tick_spc) # every 100 bins
    axmatrix.set_xticks(xticks)
    axmatrix.set_xticklabels(xticks * step_size/1000)
    axmatrix.set_yticks(range(D.shape[0]))
    axmatrix.set_yticklabels(np.array(names)[idx], minor=False, )
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()

    # Plot colorbar.
    fig.colorbar(im, cax=axcolor, orientation='horizontal')

    if save_plot:
        plt.savefig('{}_hclust.pdf'.format(outname), bbox_inches='tight')
        pd.DataFrame(D.T, columns=names).to_csv('{}_dat.csv'.format(outname))
    else:
        plt.show()

    plt.close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--use-polo', action='store_true', default=False, help='Enable optimal leaf ordering')
    parser.add_argument('-b', '--bin-size', default=5000, type=int, help='Bin size in BP [5000]')
    parser.add_argument('-t', '--tick-spacing', default=100, type=int, help='Tick spacing in bins [100]')
    parser.add_argument('input_dir', help='Input directory containing coverage files')
    parser.add_argument('output_file', help='Output filename')
    args = parser.parse_args()

    print 'Reading and binning coverage files from {}'.format(args.input_dir)
    meds, bins = get_block_medians(args.input_dir, step_size=args.bin_size)

    print 'Clustering and plotting...'
    run_hclust(args.output_file, meds, bins[0], args.bin_size, args.tick_spacing, save_plot=True, use_polo=args.use_polo)

