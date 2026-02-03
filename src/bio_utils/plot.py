"""Functions for plotting."""

from collections import defaultdict
from math import inf, log2

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.text import TextPath
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Affine2D

import bio_utils.color as bu_color
import bio_utils.sequence as bu_sequence


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
        colormap = bu_color.protein_color_map

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
    try:
        alphabet = set(alphabet)
    except TypeError:
        raise RuntimeError('Failed to coerce parameter alphabet into a set.')
    if len(alphabet) == 0:
        raise RuntimeError('Argument alphabet is empty.')
    if not all([isinstance(sym, str) and len(sym) == 1 for sym in alphabet]):
        raise RuntimeError('Not all elements of argument alphabet are single character strings.')

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
    alphabet_records = {
        'protein': (bu_sequence.protein_alphabet_strict, bu_color.protein_colormap),
        'nucleic': (bu_sequence.nucleic_alphabet_strict, bu_color.nucleic_colormap),
    }
    if alphabet in alphabet_records:
        alphabet, default_colormap = alphabet_records[alphabet]
        if colormap is None:
            colormap = default_colormap
    else:
        raise ValueError(
            f"'{alphabet}' is not a valid value for alphabet; "
            f"supported values are 'protein' and 'nucleic'"
        )
    if ax is None:
        _, ax = plt.subplots()
    if colormap is None:
        colormap = bu_color.protein_color_map
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


def plot_tree(
    tree,
    *,
    ax=None,
    linecolor=None,
    linewidth=None,
    tip_labels=True,
    tip_fontname=None,
    tip_fontsize=None,
    tip_fontcolor=None,
    tip_offset=0.005,
    support_labels=False,
    support_format_spec=None,
    support_fontname=None,
    support_fontsize=None,
    support_fontcolor=None,
    support_ha='center',
    support_va='top',
    support_hoffset=0,
    support_voffset=0,
    xmin_pad=0.01,
    xmax_pad=0.1,
    ymin_pad=0.025,
    ymax_pad=0.025,
):
    """Plot tree.

    Parameters
    ----------
    tree : TreeNode
    ax : Axes
        Axes used to draw tree. If None, a new Figure and Axes are created.
    linecolor : color (matplotlib) or dict
        If is dict, linecolor is a mapping of nodes to colors.
    linewidth : float
        Width of branches.
    tip_labels : bool
        Toggle drawing tip labels.
    tip_fontname : str
    tip_fontsize : float or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
    tip_fontcolor : color (matplotlib)
    tip_offset : float
        The offset of the tip label from the end of the branch tip in axes
        coordinates.
    support_labels : bool
        Toggle drawing support labels. Support labels are obtained from
        support attribute. If the attribute is None, it is ignored. Use
        assign_supports to extract the numerical values from the node labels.
    support_format_spec : str
        Format specification for supports using the format specification mini-
        language.
    support_fontname : str
    support_fontsize : float
    support_fontcolor : color (matplotlib)
    support_ha : {'left', 'right', 'center'}
        The horizontal alignment of the support label relative to the branch.
    support_va : {'center', 'top', 'bottom', 'baseline', 'center_baseline'}
        The vertical alignment of the support label relative to the branch.
    support_voffset : float
        The vertical offset of the support label relative to its alignment in
        axes coordinates.
    support_hoffset : float
        The horizontal offset of the support label relative to its alignment in
        axes coordinates.
    xmin_pad, xmax_pad, ymin_pad, ymax_pad : float
        Fraction of respective data ranges to pad on lower and upper bounds of
        respective axes.

    Returns
    -------
    ax : Axes
    """
    # Set options
    if ax is None:
        _, ax = plt.subplots()
    if not isinstance(linecolor, dict):
        if linecolor is None:
            linecolor = 'black'
        node2color = defaultdict(lambda x=linecolor: x)  # Use default parameter to ensure lambda captures value # fmt: skip
    else:
        node2color = linecolor
    if linewidth is None:
        linewidth = plt.rcParams['lines.linewidth']
    if tip_fontname is None:
        tip_fontname = plt.rcParams['font.family']
    if tip_fontsize is None:
        tip_fontsize = plt.rcParams['font.size']
    if tip_fontcolor is None:
        tip_fontcolor = 'black'
    if support_fontname is None:
        support_fontname = plt.rcParams['font.family']
    if support_fontsize is None:
        support_fontsize = plt.rcParams['font.size']
    if support_fontcolor is None:
        support_fontcolor = 'black'
    if support_ha == 'left':
        get_x = lambda node: xpos.get(node.parent, 0)  # In case root node
    elif support_ha == 'center':
        get_x = lambda node: (xpos.get(node.parent, 0) + xpos[node]) / 2
    elif support_ha == 'right':
        get_x = lambda node: xpos[node]
    else:
        raise ValueError(
            f"'{support_ha}' is not a valid value for support_ha; "
            f"supported values are 'left', 'right', 'center'"
        )

    # Make mapping of each node to its horizontal positions
    xpos = {}
    for node in tree.traverse(order='pre'):  # Return nodes on the way in
        length = node.length if node.length is not None else 0
        depth = xpos[node.parent] if node.parent else 0  # Checks for root node
        xpos[node] = depth + length

    # Make mapping of each node to its vertical position
    tips = list(tree.tips())
    ymax = len(tips) - 1
    ypos = {tip: ymax - i for i, tip in enumerate(reversed(tips))}
    for node in tree.traverse(order='post'):  # Return nodes on the way out
        if node.children:
            ypos[node] = (ypos[node.children[0]] + ypos[node.children[-1]]) / 2

    # Adjust axes and add labels
    xmin, xmax = min(xpos.values()), max(xpos.values())
    xrange = xmax - xmin
    ax.set_xlim(xmin - xmin_pad * xrange, xmax + xmax_pad * xrange)

    ymin, ymax = min(ypos.values()), max(ypos.values())
    yrange = ymax - ymin
    ax.set_ylim(ymax + ymax_pad * yrange, ymin - ymin_pad * yrange)  # Invert the y-axis (origin at the top) # fmt: skip

    # Plot lines and text
    lines = []
    T = ax.transLimits.inverted()  # Axes to data coordinates
    for node in tree.traverse(order='post'):  # So parent lines are drawn over children lines
        linecolor = node2color[node]
        if node.parent:  # Horizontal line of node
            x0, x1 = xpos[node.parent], xpos[node]
            y = ypos[node]
            lines.append(([(x0, y), (x1, y)], linewidth, linecolor))
        if node.children:  # Vertical line of node
            x = xpos[node]
            y0, y1 = ypos[node.children[-1]], ypos[node.children[0]]
            lines.append(([(x, y0), (x, y1)], linewidth, linecolor))
        if tip_labels and node.is_tip():  # Write tip names
            dx, _ = T.transform((tip_offset, 0)) - T.transform((0, 0))
            ax.text(
                xpos[node] + dx,
                ypos[node],
                node.name,
                verticalalignment='center',
                fontname=tip_fontname,
                fontsize=tip_fontsize,
                color=tip_fontcolor,
            )
        if (support_labels and node.support is not None and not node.is_tip()):  # Write support values if not None and not tip # fmt: skip
            if support_format_spec is None:
                support_string = str(node.support)
            else:
                support_string = f'{node.support:{support_format_spec}}'
            dx, dy = T.transform((support_hoffset, support_voffset)) - T.transform((0, 0))
            ax.text(
                get_x(node) + dx,
                ypos[node] + dy,
                support_string,
                fontname=support_fontname,
                fontsize=support_fontsize,
                color=support_fontcolor,
                horizontalalignment=support_ha,
                verticalalignment=support_va,
            )
    lc_args = {key: value for key, value in zip(['segments', 'linewidth', 'color'], zip(*lines))}
    ax.add_collection(LineCollection(**lc_args, capstyle='projecting'))

    return ax
