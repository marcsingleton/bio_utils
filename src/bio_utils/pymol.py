"""Functions for PyMOL visualizations."""

import fnmatch
import glob
from itertools import groupby

import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymol
import pymol.cmd as cmd
import pymol.cgo as cgo

import bio_utils.color as bu_color


_plldt_records = [
    ('very_low', [255, 125, 69], 50),
    ('low', [255, 219, 19], 70),
    ('high', [101, 203, 243], 90),
    ('very_high', [0, 83, 214], 100),
]

_hydrophobicity_scales = {}

_kyte_doolittle_records = [
    ('A', 1.8),
    ('C', 2.5),
    ('D', -3.5),
    ('E', -3.5),
    ('F', 2.8),
    ('G', -0.4),
    ('H', -3.2),
    ('I', 4.5),
    ('K', -3.9),
    ('L', 3.8),
    ('M', 1.9),
    ('N', -3.5),
    ('P', -1.6),
    ('Q', -3.5),
    ('R', -4.5),
    ('S', -0.8),
    ('T', -0.7),
    ('V', 4.2),
    ('W', -0.9),
    ('Y', -1.3),
]
_hydrophobicity_scales['kyte-doolittle'] = _kyte_doolittle_records

_eisenberg_records = [
    ('A', 0.62),
    ('C', 0.29),
    ('D', -0.9),
    ('E', -0.74),
    ('F', 1.19),
    ('G', 0.48),
    ('H', -0.4),
    ('I', 1.38),
    ('K', -1.5),
    ('L', 1.06),
    ('M', 0.64),
    ('N', -0.78),
    ('P', 0.12),
    ('Q', -0.85),
    ('R', -2.53),
    ('S', -0.18),
    ('T', -0.05),
    ('V', 1.08),
    ('W', 0.81),
    ('Y', 0.26),
]
_hydrophobicity_scales['eisenberg'] = _eisenberg_records


def color_plddt(selection='*'):
    """
    Color atoms with the AlphaFold-style pLDDT scheme using the stored B-factors.

    Parameters
    ----------
    selection : str
        Selection-expression or name-pattern corresponding to the atoms or objects to be colored.
    """
    b_lower = 0
    for name, rgb, b_upper in _plldt_records:
        color = '0x' + bu_color.to_hex(rgb).removeprefix('#')
        cmd.color(
            color, f'({selection}) and (b > {b_lower}) and ((b < {b_upper}) or (b = {b_upper}))'
        )
        b_lower = b_upper


cmd.extend(color_plddt)


def color_hydrophobicity(selection='*', scale='kyte-doolittle', cmap='PiYG'):
    """
    Color atoms by hydrophobicity.

    Parameters
    ----------
    selection : str
        Selection-expression or name-pattern corresponding to the atoms or objects to be colored.
    scale : str
        Name of registered hydrophobicity scale.
    cmap : Colormap or str
        Matplotlib Colormap or name of registered cmap.
    """
    if scale not in _hydrophobicity_scales:
        raise RuntimeError('Argument scale is not a registered hydrophobicity scale.')
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    elif not isinstance(mpl_colors.Colormap):
        raise RuntimeError('Argument cmap is not the name of a registered Colormap or a Colormap.')

    records = _hydrophobicity_scales[scale]
    vmin = min(record[1] for record in records)
    vmax = max(record[1] for record in records)
    norm = plt.Normalize(vmin, vmax)

    for sym, value in records:
        color = '0x' + mpl_colors.to_hex(cmap(norm(value))).removeprefix('#')
        cmd.color(color, f'({selection}) and pepseq {sym}')


cmd.extend(color_hydrophobicity)


def color_from_table(
    value_field,
    data=None,
    data_path=None,
    sep=None,
    selection='*',
    id_field='id',
    cmap='Reds',
    normalize=True,
    read_table_kwargs=None,
):
    """
    Color atoms by values in a dataframe or CSV file.

    value_field : str
        The column containing the values.
    data: DataFrame
        A dataframe containing the data loaded into memory.
    data_path : str
        Path to the CSV file.
    sep : str
        Field separator. If None, will attempt to infer from extension (.tsv or .csv). If inference
        fails, will defer to the pandas parsing engine.
    selection : str
        Selection-expression or name-pattern corresponding to the atoms or objects to be colored.
    id_field: str
        The column containing the atom id.
    cmap : Colormap or str
        Matplotlib Colormap or name of registered cmap.
    normalize : bool
        If True, scale the values with min-max transformation.
    read_table_kwargs : dict
        Additional arguments unpacked in read_table call to load data from path.
    """
    if data is not None and data_path is not None:
        raise RuntimeError('Arguments data and data_path cannot both be provided.')
    elif data is not None:
        pass  # Already assigned to data
    elif data_path is not None:
        if sep is not None:
            pass
        elif data_path.endswith('.tsv'):
            sep = '\t'
        elif data_path.endswith('.csv'):
            sep = ','

        if read_table_kwargs is None:
            read_table_kwargs = {}
        elif not isinstance(read_table_kwargs, dict):
            raise RuntimeError('Argument read_table_kwargs is not a dict.')

        data = pd.read_table(
            data_path,
            sep=sep,
            usecols=[id_field, value_field],
            **read_table_kwargs,
        )
    else:
        raise RuntimeError('Neither argument data or data_path was provided.')

    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    elif not isinstance(mpl_colors.Colormap):
        raise RuntimeError('Argument cmap is not the name of a registered Colormap or a Colormap.')

    if normalize:
        vmin, vmax = data[value_field].min(), data[value_field].max()
        values = (data[value_field] - vmin) / (vmax - vmin)
    else:
        values = data[value_field]

    for atom_id, value in zip(data[id_field], values):
        color = cmap(value)
        color_string = '0x' + mpl_colors.to_hex(color)[1:]  # Trim preceding #
        cmd.color(color_string, f'({selection}) and id {atom_id}')


cmd.extend(color_from_table)


def set_style(style_name, styles=None):
    """
    Apply a pre-defined style.

    Parameters
    ----------
    style_name : str
        Name of a pre-defined style.
    styles : dict
        Nested dict where outer keys are style names and inner dicts are setting-value pairs. If is
        None, attempts to use a STYLES variable defined in the PyMOL CLI.
    """
    if styles is None:
        if 'STYLES' in pymol.__dict__:
            styles = pymol.__dict__['STYLES']
        else:
            raise RuntimeError('Argument style not provided and STYLES not defined in PyMOL CLI.')
    elif not isinstance(styles, dict):
        raise RuntimeError('Argument styles is not a dict.')

    try:
        style = styles[style_name]
    except KeyError:
        raise RuntimeError(f'Style {style_name} not recognized')

    for setting, value in style.items():
        cmd.set(setting, value)


cmd.extend(set_style)


def load_color_scheme(name, as_cycle=False):
    """
    Add color scheme from bio_utils.color module into color namespace.

    Parameters
    ----------
    name : str
        Name of color scheme. See bio_utils.color for details.
    as_cycle : bool
        If True, names colors with as 'C{index}'. Otherwise named as '{name}_{index}'.
    """
    # Some bookkeeping for internal clarity
    scheme_name = name
    del name

    if scheme_name not in bu_color.color_schemes:
        raise RuntimeError(f'Color scheme {scheme_name} not found in defined color schemes.')
    if as_cycle:
        prefix = 'C'
    else:
        prefix = f'{scheme_name}_'

    color_scheme = bu_color.color_schemes[scheme_name]
    for i, color in enumerate(color_scheme):
        cmd.set_color(f'{prefix}{i}', mpl_colors.to_rgb(color))


cmd.extend(load_color_scheme)


def load_glob(pattern, recursive=False, load_fn=None, name_fn=None):
    """
    Load multiple files matching a glob pattern.

    Parameters
    ----------
    pattern : str
        Glob pattern to match files (e.g., '*.pdb').
    recursive : bool, optional
        If True, the pattern '**' will match files and directories recursively.
    load_fn : function
        Function to load file at path. Accepts a path string as its first argument and, optionally,
        a name for the loaded object as its second argument.
    name_fn : function
        Function to generate names for loaded objects. Accepts a path string as a single argument
        and returns a string.
    """
    if load_fn is None:
        load_fn = cmd.load
    paths = glob.glob(pattern, recursive=recursive)
    if name_fn is None:
        for path in paths:
            load_fn(path)
    else:
        for path in paths:
            name = name_fn(path)
            load_fn(path, name)


cmd.extend(load_glob)


def load_cgo_arrow(
    origin,
    orient,
    name,
    color,
    tail_length=5,
    tail_width=0.5,
    head_length=1.5,
    head_width=1.5,
    mode='tail',
    pseudoatom=False,
):
    """
    Load 3D arrow as object.

    Parameters
    ----------
    origin : list of three floats
        Vector of position of arrow.
    orient : list of three floats
        Vector of direction of arrow head. Does not need to be a unit vector.
    name : str
        Name of arrow object.
    color : str or list of three floats
        Registered color name or RGB values in range [0.0, 1.0] or [0, 255].
    tail_length : float
        Length of tail in A.
    tail_width : float
        Width of tail in A.
    head_length : float
        Length of head in A.
    head_width : float
        Width of head in A.
    mode : tail | tail-center | center
        Sets alignment of arrow with origin. If tail, arrow tail starts at origin. tail-center and
        center align those points with the origin.
    pseudoatom : bool
        If True, add pseudoatom at arrow tip as new object with provided name and '_pk' suffix.

    PyMOL seemingly cannot accommodate CGO objects and pseudoatoms under a single name, so they must
    be separate objects.
    """
    if isinstance(color, str):
        index = cmd.get_color_index(color)
        if index == -1:
            raise RuntimeError(f'Color {color} is not defined')
        r, g, b = cmd.get_color_tuple(index)
    else:
        r, g, b = color
    origin = np.asarray(origin)
    orient = np.asarray(orient) / np.linalg.vector_norm(orient)

    accepted_modes = ['tail', 'tail-center', 'center']
    if mode == 'tail':
        pass
    elif mode == 'tail-center':
        length = tail_length
        origin -= length / 2 * orient
    elif mode == 'center':
        length = tail_length + head_length
        origin -= length / 2 * orient
    else:
        raise ValueError(f'Mode {mode} not recognized. Must be in {accepted_modes}.')

    # fmt: off
    t0, t1 = origin, origin + tail_length * orient
    tail = [
        cgo.CYLINDER,
        t0[0], t0[1], t0[2],
        t1[0], t1[1], t1[2],
        tail_width,
        r, g, b,
        r, g, b,
    ]
    h0, h1 = t1, t1 + head_length * orient
    head = [cgo.CONE,
            h0[0], h0[1], h0[2],
            h1[0], h1[1], h1[2],
            head_width, 0,
            r, g, b,
            r, g, b,
            1, 1,
        ]
    arrow = tail + head
    # fmt: on

    view = cmd.get_view()
    cmd.load_cgo(arrow, name)
    if pseudoatom:
        pos = origin + (tail_length + head_length) * orient
        cmd.pseudoatom(name + '_pk', pos=tuple(pos))  # Conversion needed for correct serialization
    cmd.set_view(view)


cmd.extend(load_cgo_arrow)


def load_vhh_arrow(
    name,
    selection='*',
    color=None,
    state=-1,
    min_strand_len=5,
    N_term_len=2,
    mode='strand-bounds',
    arrow_kwargs=None,
):
    """
    Create an arrow beginning at the beta sandwich and pointing along its major axis.

    Parameters
    ----------
    name : str
        Name of arrow object.
    selection : str
        Selection-expression or name-pattern corresponding to the VHH atoms or object.
    color : None or str or list of three floats
        If None, infer from first carbon atom in selection. Otherwise, a registered color name or
        RGB values in range [0.0, 1.0] or [0, 255].
    state : int
        State containing secondary structure annotations.
    min_strand_len : int
        Minimum number of residues in a beta strand.
    N_term_len : int
        Number of residues at N-terminus used for identifying CDR orientation.
    mode : 'strand-bounds' | 'pca'
        Method for determining arrow origin and orientation.
    arrow_kwargs : dict
        Additional arguments unpacked in load_cgo_arrow call.
    """
    if min_strand_len < 1:
        raise ValueError('min_strand_len must be greater than 1')
    if N_term_len < 2:
        raise ValueError('N_term_len must be greater than 2')
    accepted_modes = ['strand-bounds', 'pca']
    if mode not in accepted_modes:
        raise ValueError(f'Mode {mode} not recognized. Must be in {accepted_modes}.')
    if arrow_kwargs is None:
        arrow_kwargs = {}
    if 'mode' not in arrow_kwargs:
        arrow_kwargs['mode'] = 'center'

    # Set color
    if color is None:
        indices = []
        cmd.iterate(f'first (({selection}) and elem C)', lambda x: indices.append(x.color))
        if not indices:
            raise RuntimeError('No carbon atoms in selection to infer color.')
        color = cmd.get_color_tuple(indices[0])

    # Gather CA atoms (with state assignments)
    atoms = []
    cmd.iterate_state(
        state,
        f'({selection}) and name CA',
        lambda atom: atoms.append((atom.ss, (atom.x, atom.y, atom.z))),
    )
    if len(atoms) == 0:
        raise RuntimeError('No atoms in selection.')

    # Check for at least one beta strand
    states = set()
    for ss, _ in atoms:
        states.add(ss)
    if 'S' not in states:
        raise RuntimeError('No beta strands in nanbody')

    # Calculate origin and orient vectors
    grouped = groupby(atoms, key=lambda x: x[0])
    if mode == 'strand-bounds':
        beta_coords = []
        boundary_pairs = []
        for key, group in grouped:
            group = [x[1] for x in group]
            if key == 'S' and len(group) >= min_strand_len:
                beta_coords.extend(group)
                boundary_pairs.append((group[0], group[-1]))
        if len(beta_coords) == 0:
            raise RuntimeError('No beta strands passed length filter.')

        # Calculate strand unit vectors
        vecs = []
        for start, stop in boundary_pairs:
            start, stop = np.array(start), np.array(stop)
            vec = stop - start
            norm = np.linalg.norm(vec)
            if np.isclose(norm, 0):
                raise RuntimeError('Length 0 vector detected. Check for beta strands of length 1.')
            vec /= np.linalg.norm(vec)
            vecs.append(vec)

        # Point vectors in same direction as first
        if len(vecs) > 1:
            ref = vecs[0]
            new_vecs = [ref]
            for vec in vecs[1:]:
                if np.dot(ref, vec) < 0:
                    new_vecs.append(-vec)
            vecs = new_vecs

        # Average coordinates
        orient = np.stack(vecs).mean(axis=0)
        origin = np.array(beta_coords).mean(axis=0)
    elif mode == 'pca':
        beta_coords = []
        for key, group in grouped:
            group = [x[1] for x in group]
            if key == 'S' and len(group) >= min_strand_len:
                beta_coords.extend(group)
        if len(beta_coords) == 0:
            raise RuntimeError('No beta strands passed length filter.')

        X = np.array(beta_coords)
        mu = X.mean(axis=0)
        X -= mu
        cov = X.T @ X
        eig_vals, eig_vecs = np.linalg.eig(cov)
        eig_argmax = np.argmax(eig_vals)

        origin = mu
        orient = eig_vecs[eig_argmax]
    else:
        raise RuntimeError(
            f'Mode {mode} not recognized. '
            'This error should be unreachable--ensure guard clause includes all recognized modes.'
        )

    # Check for proper orientation
    if N_term_len > len(atoms):
        N_term_len = len(atoms)
        raise RuntimeWarning(
            'N_term_len is greater than selection length. Setting to selection length.'
        )
    N_term_coord_0 = np.array(atoms[0][1])
    N_term_coord_1 = np.array(atoms[N_term_len - 1][1])
    N_term_vec = N_term_coord_1 - N_term_coord_0
    if np.dot(N_term_vec, orient) > 0:
        orient *= -1  # CDRs point in opposite direction as N-terminus

    load_cgo_arrow(origin, orient, name, color, **arrow_kwargs)


cmd.extend(load_vhh_arrow)


def load_as_vhh_arrow(
    path,
    name=None,
    selection='*',
    color=None,
    arrow_kwargs=None,
    vhh_kwargs=None,
):
    """
    Load a VHH structure as an arrow.

    Parameters
    ----------
    path : str
        Path to structure.
    name : str
        Name of arrow object.
    selection : str
        Selection-expression or name-pattern corresponding to the VHH atoms or object.
    color : None or str or list of three floats
        If None, infer from first carbon atom in selection. Otherwise, a registered color name or
        RGB values in range [0.0, 1.0] or [0, 255].
    arrow_kwargs : dict
        Additional arguments passed to load_vhh_arrow call.
    vhh_kwargs : dict
        Additional arguments unpacked in load_vhh_arrow call.
    """
    # Some bookkeeping for internal clarity
    arrow_name = name
    del name

    if vhh_kwargs is None:
        vhh_kwargs = {}
    if arrow_kwargs is None:
        arrow_kwargs = {}

    cmd.load(path)
    model_name = cmd.get_object_list()[-1]

    # Set name and color
    if arrow_name is None:
        arrow_name = model_name
    if arrow_name == model_name:
        tmp_name = cmd.get_unused_name()
        cmd.set_name(model_name, tmp_name)
        model_name = tmp_name

    load_vhh_arrow(
        arrow_name,
        color,
        f'({selection}) and %{model_name}',
        arrow_kwargs=arrow_kwargs,
        **vhh_kwargs,
    )
    cmd.delete(f'%{model_name}')


cmd.extend(load_as_vhh_arrow)


def apply(fn, pattern='*', selection='*'):
    """Apply a function to a selection in a list of objects.

    Parameters
    ----------
    fn : function
        Function to apply. Accepts the object name and selection as the first and second parameters,
        respectively.
    pattern : str
        Name-pattern corresponding to objects.
    selection : str
        Selection-expression to be applied to each object.
    """
    for name in cmd.get_object_list():
        if not fnmatch.fnmatch(name, pattern):
            continue
        fn(name, selection)


cmd.extend(apply)
