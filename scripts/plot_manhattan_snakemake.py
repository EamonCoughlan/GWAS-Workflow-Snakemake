#this is the version adapted by Eamon for Snakemake, from Konrad's original 'plot_manhatten.py'
import sys
import pandas as pd
from numpy import *

from matplotlib.pyplot import *
import seaborn as sb


def genomic_pos(pos, col_chr, col_bp):
    l = pos.groupby(col_chr)[col_bp].max()
    offset = l.cumsum() - l
    return pos[col_bp] + list(offset.loc[pos[col_chr]])


def chromosome_colors(N, name='Alternating'):
    colors = []
    if name == 'Alternating':
        colors = sb.color_palette("Paired", 6)
        colors = [colors[1], colors[5]]
        colors = (colors * int(ceil(N / 2)))[:N]
    elif name == 'Alternating-bw':
        colors = [(0.4, 0.4, 0.4), (.0, .0, .0)]
        colors = (colors * int(ceil(N / 2)))[:N]
    else:
        colors = sb.color_palette(name, N)
    assert (len(colors) == N)
    return colors


def guess_col_chr(cols):
    what = list(filter(lambda x: x in cols, ['Chr', 'CHR', 'chrom', 'Chrom', 'CHROM', 'chromosome', '#CHROM']))
    if 1 == len(what):
        return what[0]
    else:
        raise ValueError(f"Could not identify chromosome column. Candidates are {what}.")


def guess_col_bp(cols):
    what = list(filter(lambda x: x in cols, ['BP', 'Bp', 'Pos', 'pos', 'position']))
    if 1 == len(what):
        return what[0]
    else:
        raise ValueError(f"Could not identify base-pair column. Candidates are {what}.")


def manhattan(values, meta, col_chr=None, col_bp=None, ylabel=None, ax=None, ylim=None, colors=None, fontsize=15,
              annotate=None):
    meta = meta.loc[values.index].copy()
    if col_chr is None:
        col_chr = guess_col_chr(meta.columns)
    if col_bp is None:
        col_bp = guess_col_bp(meta.columns)

    meta['genomic_position'] = genomic_pos(meta, col_chr, col_bp)
    data = pd.concat([meta, values.to_frame('value')], axis=1)

    Nchr = len(data[col_chr].unique())
    if colors is None:
        colors = chromosome_colors(Nchr)
    elif type(colors) is str:
        colors = chromosome_colors(Nchr, name=colors)
    elif len(colors) == N:
        pass
    else:
        raise ValueError(f"Argument colors with value {colors} not understood.")

    # compute ticks:
    ticks = {}
    for c, g in data.groupby(col_chr):
        bounds = [g.genomic_position.min(), g.genomic_position.max()]
        ticks[c] = bounds[0] + (bounds[1] - bounds[0]) / 2

    if ax is None:
        fig, ax = subplots(figsize=(15, 8))

    which = data.index if ylim is None else (ylim[0] <= data.value) & (data.value <= ylim[1])
    for c, g in data.loc[which].groupby(col_chr):
        ax.scatter(g.genomic_position, g.value, c=np.array([colors[c - 1]]), marker='.', rasterized=True, )
    if ylim is not None:
        which = ylim[0] > data.value
        for c, g in data.loc[which].groupby(col_chr):
            ax.scatter(g.genomic_position, [ylim[0] - 2] * g.shape[0], c=np.array([colors[c - 1]]), marker='|',
                       rasterized=True, )
        which = ylim[1] < data.value
        for c, g in data.loc[which].groupby(col_chr):
            ax.scatter(g.genomic_position, [ylim[1] + 2] * g.shape[0], c=np.array([colors[c - 1]]), marker='|',
                       rasterized=True, )
    if annotate is not None:
        l = data.groupby(col_chr)[col_bp].max()
        offset = l.cumsum() - l
        y_top = ax.get_ylim()[1] + 1
        for name in annotate:
            reg = annotate[name]
            l = offset.loc[reg[0]] + reg[1]
            u = offset.loc[reg[0]] + reg[2]
            ax.plot([l, l, u, u], [y_top - 1, y_top, y_top, y_top - 1], '-', color='black')
            ax.annotate(xy=[l + (u - l) / 2., y_top + 0.5], xytext=[0, 1], textcoords='offset points', xycoords='data',
                        s=str(name), annotation_clip=False, ha='center')
    ax.set_xticks(list(ticks.values()))
    ax.set_xticklabels(list(ticks.keys()), fontsize=fontsize)
    ax.set_xlabel('Genomic Position', fontsize=fontsize)
    ax.set_title(snakemake.wildcards.group)
    if not ylabel is None:
        ax.set_ylabel(ylabel, fontsize=fontsize)


def main(infile):
    files = []
    with open(infile, 'r') as list_of_files:
        for file in list_of_files :
            files.append(file.strip('\n'))
    dfs = [pd.read_csv(f, sep='\s+') for f in files]
    df = pd.concat(dfs).reset_index(drop=True).sort_values(['CHROM'])

    fig, ax = subplots(1, figsize=(15, 6))
    manhattan(df['LOG10P'], df, col_chr='CHROM', col_bp='GENPOS', ax=ax, ylabel='$-\log_{10} P$')
    sb.despine(offset=19, trim=True)
    fig.savefig(snakemake.output[0], format='png', bbox_inches='tight')

main(snakemake.input[0])
