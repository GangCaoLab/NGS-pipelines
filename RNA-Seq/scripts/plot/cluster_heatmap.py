import click
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


def re_arrange_first2col(X, y):
    labels = np.unique(y)
    parts = []
    for i in labels:
        subpart = X[y==i]
        parts.append(subpart)
    #parts.sort(key=lambda p: p[:, 0].mean())
    res = np.concatenate(parts, axis=0)
    return res


def kmeans_arrange(m, n_clusters):
    X = m.as_matrix()
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X)
    y = kmeans.labels_
    X = re_arrange_first2col(X, y)
    m = pd.DataFrame(X, columns=m.columns)
    return m


@click.command()
@click.argument("file_in")
@click.argument("fig_out")
@click.option("--cluster-type", "-t",
    type=click.Choice(['hierachical', 'kmeans']),
    default='hierachical',
    help="The algorithm for cluster genes.")
@click.option("--k-num", "-k",
    type=int,
    default=12,
    help="How many clusters(k-number) for K-Means algorithm")
@click.option("--sample",
    type=int,
    default=-1,
    help="Down sampling the input table, "
         "e.g. --sample 1000, will select 1000 genes randomly to plot.")
@click.option("--col-c/--no-col-c", "colc",
    default=False,
    help="cluster columns or not. default False")
@click.option("--keep-gene-name/--no-gene-name", "keep_name",
    default=False,
    help="Show gene's name or not.")
@click.option("--zscore/--no-zscore", "zscore",
    default=True,
    help="use zscore replace real value or not.")
@click.option("--cmap",
    default="RdYlBu_r",
    help="The color map of heatmap, default RdYlBu_r")
@click.option("--dpi",
    default=600,
    help="The dpi of output figure.")
def cluster_heatmap(
        file_in, fig_out,
        cluster_type, k_num,
        sample,
        colc, keep_name, zscore,
        cmap, dpi):
    """
    plot cluster heatmap, support K-Means and Hierachical cluster algorithm.

    \b
    NOTE:
        1. input file, must be tab delimited table
    """
    # read table
    m = pd.read_csv(file_in, sep='\t', index_col=0)
    m.index.name = 'gene'

    # down sampling genes
    if sample != -1:
        m = m.sample(n=sample)

    # rearange all rows if use K-Means
    if cluster_type == 'kmeans':
        m = kmeans_arrange(m, k_num)
        row_c = False
    elif cluster_type == 'hierachical':
        row_c = True

    # draw cluster heatmap
    z = 0 if zscore else None
    ch_m = sns.clustermap(m, cmap=cmap, row_cluster=row_c, col_cluster=colc, z_score=z)

    if not keep_name: # remove gene's name
        ch_m.ax_heatmap.set_yticks([])
        ch_m.ax_heatmap.set_yticklabels([])

    # draw divider line between columns
    for i in range(1, m.shape[1]):
        ch_m.ax_heatmap.axvline(i, c='w', linewidth=0.2)

    # save figure
    plt.savefig(fig_out, dpi=dpi)


if __name__ == "__main__":
    cluster_heatmap()
