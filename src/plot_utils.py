"""Functions for plotting."""

from math import inf, log2

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.text import TextPath
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Affine2D

tableau_10 = [
    '#4E79A7',  # blue
    '#F28E2B',  # orange
    '#E15759',  # red
    '#76B7B2',  # cyan
    '#59A14F',  # green
    '#EDC948',  # yellow
    '#B07AA1',  # purple
    '#FF9DA7',  # pink
    '#9C755F',  # brown
    '#BAB0AC',  # grey
]
tableau_20 = [
    '#4E79A7',  # blue
    '#F28E2B',  # orange
    '#E15759',  # red
    '#499894',  # dark_cyan
    '#59A14F',  # green
    '#B6992D',  # dark_yellow
    '#B07AA1',  # purple
    '#D37295',  # maroon
    '#9D7660',  # dark_brown
    '#79706E',  # dark_grey
    '#A0CBE8',  # light_blue
    '#FFBE7D',  # light_orange
    '#FF9D9A',  # light_red
    '#86BCB6',  # light_cyan
    '#8CD17D',  # light_green
    '#F1CE63',  # light_yellow
    '#D4A6C8',  # light_purple
    '#FABFD2',  # light_maroon
    '#D7B5A6',  # light_brown
    '#BAB0AC',  # grey
]

protein_colormap = {
    'A': '#6DD7A1',
    'I': '#55C08C',
    'L': '#55C08C',
    'V': '#55C08C',
    'M': '#55C08C',
    'F': '#B897EC',
    'Y': '#B897EC',
    'W': '#A180D2',
    'S': '#FFBE74',
    'T': '#FFBE74',
    'N': '#77EAF4',
    'Q': '#77EAF4',
    'D': '#EE8485',
    'E': '#EE8485',
    'H': '#96C4FF',
    'K': '#7FADEA',
    'R': '#7FADEA',
    'C': '#FAED70',
    'G': '#E2DEDD',
    'P': '#FFB1F1',
    'X': '#93908F',
    '-': '#FFFFFF',
    '.': '#3F3F3F',
}
nucleic_colormap = {
    'A': '#6DD7A1',
    'T': '#EE8485',
    'G': '#E2DEDD',
    'C': '#7FADEA',
    'N': '#93908F',
}


def _get_ratio(rows, cols, hspace, data_height, data_hspace, max_cols):
    block_num, _, block_cols = _get_block_dims(rows, cols, max_cols)
    block_height = block_num * (1 + data_hspace + data_height) * rows
    hspace_height = (block_num - 1) * hspace * rows
    height = block_height + hspace_height
    width = block_cols
    ratio = width / height
    return ratio


def _get_block_dims(rows, cols, max_cols):
    if max_cols >= cols:
        block_num = 1
        block_cols = cols
    else:
        q, r = divmod(cols, max_cols)
        block_num = q + int(r > 0)
        block_cols = max_cols
    block_rows = rows
    return block_num, block_rows, block_cols


def _get_best_max_cols(rows, cols, width, height, hspace, data_height, data_hspace):
    target_ratio = width / height
    best_delta = inf
    best_max_cols = None
    for max_cols in range(cols, 0, -1):
        ratio = _get_ratio(rows, cols, hspace, data_height, data_hspace, max_cols)
        delta = abs(target_ratio - ratio)
        if delta < best_delta:
            best_delta = delta
            best_max_cols = max_cols
        elif delta >= best_delta:
            return best_max_cols


def plot_alignment(
    alignment,
    fig=None,
    hspace=0.5,
    data_axs=False,
    data_height=0.5,
    data_hspace=0.1,
    start=1,
    max_cols=None,
    show_labels=True,
    label_kwargs=None,
    show_syms=True,
    sym_kwargs=None,
    colormap=None,
):
    """
    Plot an alignment.

    Parameters
    ----------
    alignment : iterable of tuples of (label, seq)
    fig : Figure
        Figure to draw alignment on. If None, will create a new Figure.
    hspace : float
        Space between blocks of MSA and data Axes in units of MSA Axes height.
    data_axs : bool
        If True, draws a data Axes for each block.
    data_height : float
        Height of data Axes in units of MSA Axes height.
    data_hspace : float
        Space between MSA and data Axes in a block in units of MSA Axes height.
    start : int
        The starting index for the x-axis.
    max_cols : int
        Maximum number of columns to display in a block. If None, will be calculated based on figure
        size.
    show_labels : bool
        If True, displays the labels on the y-axis.
    label_kwargs : dict
        Additional keyword arguments to pass to the label text rendering.
    show_syms : bool
        If True, displays the symbols in the MSA.
    sym_kwargs : dict
        Additional keyword arguments to pass to the symbol text rendering.
    colormap : dict
        A dictionary mapping symbols to their corresponding color codes. If None, defaults to
        `protein_colormap`.

    Returns
    --------
        fig : Figure
        axs : list of Axes
    """
    if fig is None:
        fig = plt.figure()
    params = fig.subplotpars
    width = fig.bbox.width * (params.right - params.left)
    height = fig.bbox.height * (params.top - params.bottom)
    if label_kwargs is None:
        label_kwargs = {}
    if sym_kwargs is None:
        sym_kwargs = {'fontname': 'monospace'}
    if colormap is None:
        colormap = protein_colormap

    if not data_axs:
        data_height = 0
        data_hspace = 0

    alignment = list(alignment)
    labels = [labels for labels, _ in alignment]
    seqs = [seq for _, seq in alignment]
    seqlens = [len(seq) for seq in seqs]

    ROWS = len(seqs)
    COLS = max(seqlens)

    if max_cols is None:
        max_cols = _get_best_max_cols(ROWS, COLS, width, height, hspace, data_height, data_hspace)

    BLOCK_NUM, BLOCK_ROWS, BLOCK_COLS = _get_block_dims(ROWS, COLS, max_cols)

    height_ratios = [1, data_hspace + data_height]
    for block_index in range(BLOCK_NUM - 1):
        height_ratios.extend([hspace, 1, data_hspace + data_height])
    gs = fig.add_gridspec(len(height_ratios), 1, height_ratios=height_ratios, hspace=0)

    axs = []
    for block_index in range(BLOCK_NUM):
        msa_ax = fig.add_subplot(gs[3 * block_index])

        patches = []
        for i in range(BLOCK_ROWS):
            y = BLOCK_ROWS - i - 1
            for j in range(BLOCK_COLS):
                j = block_index * BLOCK_COLS + j
                if j >= len(seqs[i]):
                    break
                x = start + j

                sym = seqs[i][j]
                xy = (x - 0.5, y - 0.5)
                patch = Rectangle(xy, 1, 1, facecolor=colormap[sym], edgecolor='none')
                patches.append(patch)
                if show_syms:
                    msa_ax.text(x, y, sym, va='center_baseline', ha='center', **sym_kwargs)

        patches = PatchCollection(patches, match_original=True)
        msa_ax.add_collection(patches)

        xmin = start + block_index * BLOCK_COLS
        xmax = start + (block_index + 1) * BLOCK_COLS

        msa_ax.set_xlim(xmin - 0.5, xmax - 0.5)
        msa_ax.set_ylim(-0.5, BLOCK_ROWS - 0.5)
        msa_ax.set_aspect('equal')
        msa_ax.xaxis.set_major_locator(MaxNLocator(nbins='auto', steps=[1, 2, 4, 5, 10]))
        if show_labels:
            msa_ax.set_yticks(range(BLOCK_ROWS), labels[::-1], **label_kwargs)
        else:
            msa_ax.tick_params(labelleft=False)
        msa_ax.tick_params(left=False)

        msa_ax.spines[:].set_visible(False)
        axs.append(msa_ax)

        if data_axs:
            data_ax = msa_ax.inset_axes((0, -(data_hspace + data_height), 1, data_height))
            data_ax.set_xlim(xmin - 0.5, xmax - 0.5)
            data_ax.xaxis.set_major_locator(MaxNLocator(nbins='auto', steps=[1, 2, 4, 5, 10]))
            msa_ax.tick_params(bottom=False, labelbottom=False)
            label_ax = data_ax
            axs.append(data_ax)
        else:
            xticks = msa_ax.get_xticks()
            label_ax = msa_ax

        if block_index >= BLOCK_NUM - 1:
            vmin = xmin
            vmax = max(seqlen + start for seqlen in seqlens) - 1
            xticks = [v for v in label_ax.get_xticks() if vmin <= v <= vmax]
            label_ax.set_xticks(xticks)

    return fig, axs


def _get_profile(alignment, alphabet, norm=True):
    """
    Calculate a profile.

    Sequences do not necessarily need to be the same length, but ragged ends are taken to be filled
    with "out of alphabet" symbols.

    Parameters
    ----------
    alignment : iterable of tuples of (label, seq)
    alphabet : iterable of single character strings
    norm: bool
        Normalize outputs to fractions.

    Returns
    -------
    p_array : list of dict[str, numeric]
        Counts or fractions of symbols at each position in alignment.
    fs : list of numeric
        Counts or fractions of "in alphabet" symbols at each position in alignment.
    """
    if alphabet == 'protein':
        # fmt: off
        alphabet = {
            'A', 'C', 'D', 'E',
            'F', 'G', 'H', 'I',
            'K', 'L', 'M', 'N',
            'P', 'Q', 'R', 'S',
            'T', 'V', 'W', 'Y',
        }
        # fmt: on
    elif alphabet == 'nucleic':
        alphabet = {'A', 'T', 'G', 'C'}
    else:
        try:
            alphabet = set(alphabet)
        except TypeError:
            raise RuntimeError('Failed to coerce parameter alphabet into a set.')
        if len(alphabet) == 0:
            raise RuntimeError('Argument alphabet is empty.')
        if not all([isinstance(str, sym) and len(sym) == 0 for sym in alphabet]):
            raise RuntimeError(
                'Not all elements of argument alphabet are single character strings.'
            )

    alignment = list(alignment)
    seqs = [seq for _, seq in alignment]
    seqlens = [len(seq) for seq in seqs]

    ROWS = len(seqs)
    COLS = max(seqlens)

    count_array = [{sym: 0 for sym in alphabet} for _ in range(COLS)]
    for seq in seqs:
        for j, sym in enumerate(seq):
            if sym in alphabet:
                count_array[j][sym] += 1

    p_array = []
    fs = []
    for counts in count_array:
        N = sum(counts.values())
        if norm:
            ps = [(sym, count / N) for sym, count in counts.items()]
            f = N / ROWS
        else:
            ps = [(sym, count) for sym, count in counts.items()]
            f = N
        p_array.append(ps)
        fs.append(f)

    return p_array, fs


def plot_profile(alignment, alphabet, ax=None, start=1, width=1, colormap=None, fontprop=None):
    """
    Plot a profile from an alignment.

    The code for scaling the symbols is adapted from the Biotite package.

    Parameters
    ----------
    alignment : iterable of tuples of (label, seq)
    alphabet : str
        Name of alphabet. 'protein' and 'nucleic' options are supported.
    ax : Axes
        Axes to draw profile on. If None, will create a New Axes.
    start : int
        The starting index for the x-axis.
    width : float
        Width of symbols.
    colormap : dict
        A dictionary mapping symbols to their corresponding color codes. If None, defaults to
        `protein_colormap`.
    fontprop : dict
        Additional parameters passed to `matplotlib.text.TextPath` for drawing symbols.

    Returns
    -------
    ax : Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    if colormap is None:
        colormap = protein_colormap
    if fontprop is None:
        {'family': ['monospace']}

    # Calculate profile
    p_array, fs = _get_profile(alignment, alphabet)
    y_array = []
    for ps, f in zip(p_array, fs):
        h = 0
        for _, p in ps:
            if p > 0:
                h -= p * log2(p)

        r = log2(len(ps)) - h

        ys = []
        for sym, p in ps:
            y = f * r * p
            ys.append((sym, y))
        y_array.append(ys)

    # Make plot
    ax.set_xlim(start - width / 2, start + len(y_array) - 1 + width / 2)
    ax.set_ylim(0, log2(len(y_array[0])))
    ax.xaxis.set_major_locator(MaxNLocator(nbins='auto', steps=[1, 2, 4, 5, 10]))
    for i, ys in enumerate(y_array):
        x = i + start
        y0 = 0
        for sym, y in sorted(ys, key=lambda x: x[1]):
            path = TextPath((0, 0), sym, size=1e-3, prop=fontprop)
            bbox = path.get_extents()
            path = (
                Affine2D()
                .translate(-bbox.x0, -bbox.y0)
                .scale(width / bbox.width, y / bbox.height)
                .translate(x - width / 2, y0)
                .transform_path(path)
            )
            patch = PathPatch(path, facecolor=colormap[sym], edgecolor='none')
            ax.add_patch(patch)
            y0 += y

    return ax
