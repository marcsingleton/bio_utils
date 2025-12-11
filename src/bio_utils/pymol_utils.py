"""Functions for PyMOL visualizations."""

import glob
from itertools import groupby

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymol.cmd as cmd
import pymol.cgo as cgo


_plldt_records = [
    ('plddt_very_low', [255, 125, 69], 50),
    ('plddt_low', [255, 219, 19], 70),
    ('plddt_high', [101, 203, 243], 90),
    ('plddt_very_high', [0, 83, 214], 100),
]


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
        cmd.set_color(name, rgb)
        cmd.color(name, f'({selection}) and (b > {b_lower}) and (b < {b_upper})')
        b_lower = b_upper


cmd.extend(color_plddt)


def color_from_table(
    value_field,
    data=None,
    data_path=None,
    sep=None,
    selection='*',
    id_field='id',
    cmap='Reds',
    normalize=True,
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
    cmap : str
        Name of cmap.
    normalize : bool
        If True, scale the values with min-max transformation.
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
        data = pd.read_table(data_path, sep=sep, usecols=[id_field, value_field])
    else:
        raise RuntimeError('Neither argument data or data_path was provided.')

    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    elif not isinstance(colors.Colormap):
        raise RuntimeError('Argument cmap is not the name of a registered Colormap or a Colormap.')

    if normalize:
        vmin, vmax = data[value_field].min(), data[value_field].max()
        values = (data[value_field] - vmin) / (vmax - vmin)
    else:
        values = data[value_field]

    for atom_id, value in zip(data[id_field], values):
        color = cmap(value)
        color_string = '0x' + colors.to_hex(color)[1:]  # Trim preceding #
        cmd.color(color_string, f'id {atom_id} and ({selection})')


cmd.extend(color_from_table)


def set_style(style_name):
    """
    Apply a pre-defined style.

    Parameters
    ----------
    style_name : str
        Name of a pre-defined style. See source for details.
    """
    if style_name == 'default':
        cmd.reinitialize('settings')
    elif style_name == 'cartoon':
        cmd.reinitialize('settings')
        cmd.bg_color('white')
        cmd.set('ray_trace_mode', 1)  # 1 renders with outlines
        cmd.set('ray_shadows', 'off')  # Turns off ray trace shadows
        cmd.set('ray_trace_disco_factor', 0.2)  # Fewer lines on surfaces
        cmd.set('antialias', 2)
        cmd.set('ambient', 0.5)  # Ambient light
        cmd.set('direct', 0.2)  # Camera light
        cmd.set('spec_direct', 0)  # Turn off camera specular reflections
        cmd.set('light_count', 1)  # Use only ambient and camera light
        cmd.set('cartoon_discrete_colors', 'on')  # Turn off color blending at boundaries
    elif style_name == 'spheres':
        cmd.reinitialize('settings')
        cmd.bg_color('white')
        cmd.set('antialias', 2)
        cmd.set('ambient', 0.2)  # Ambient light
        cmd.set('spec_direct', 0)  # Turn off camera specular reflections
        cmd.set('light_count', 1)  # Use only ambient and camera light
    elif style_name == 'surface':
        cmd.reinitialize('settings')
        cmd.bg_color('white')
        cmd.set('ray_trace_mode', 1)  # 1 renders with outlines
        cmd.set('ray_shadows', 'off')  # Turns off ray trace shadows
        cmd.set('ray_trace_disco_factor', 1)  # No lines on surfaces
        cmd.set('antialias', 2)
        cmd.set('ambient', 0.65)  # Ambient light
        cmd.set('direct', 0)  # Camera light
        cmd.set('spec_direct', 0)  # Turn off camera specular reflections
        cmd.set('light_count', 1)  # Use only ambient and camera light
        cmd.set('surface_quality', 1)
    elif style_name == 'high_contrast':
        cmd.bg_color('white')
        cmd.set('ray_trace_mode', 1)  # 1 renders with outlines
        cmd.set('ray_shadows', 'off')  # Turns off ray trace shadows
        cmd.set('ray_trace_disco_factor', 1)  # No lines on surfaces
        cmd.set('antialias', 2)
        cmd.set('ambient', 1)  # Ambient light
        cmd.set('direct', 0)  # Camera light
        cmd.set('spec_direct', 0)  # Turn off camera specular reflections
        cmd.set('light_count', 1)  # Use only ambient and camera light
        cmd.set('cartoon_discrete_colors', 'on')  # Turn off color blending at boundaries
        cmd.set('surface_quality', 1)
        cmd.set('depth_cue', 'off')
    else:
        raise RuntimeError(f'Style {style_name} not recognized')


cmd.extend(set_style)


def add_palette(palette_name):
    """
    Add colors from a pre-defined palette into color namespace.

    Parameters
    ----------
    palette_name : str
        Name of pre-defined palette. See source for details.
    """
    if palette_name == 'tableau_10':
        cmd.set_color('C0', [78, 121, 167])
        cmd.set_color('C1', [242, 142, 43])
        cmd.set_color('C2', [225, 87, 89])
        cmd.set_color('C3', [118, 183, 178])
        cmd.set_color('C4', [89, 161, 79])
        cmd.set_color('C5', [237, 201, 72])
        cmd.set_color('C6', [176, 122, 161])
        cmd.set_color('C7', [255, 157, 167])
        cmd.set_color('C8', [156, 117, 95])
        cmd.set_color('C9', [186, 176, 172])
    elif palette_name == 'tableau_20':
        cmd.set_color('C0-A', [78, 121, 167])
        cmd.set_color('C0-B', [160, 203, 232])
        cmd.set_color('C1-A', [242, 142, 43])
        cmd.set_color('C1-B', [255, 190, 125])
        cmd.set_color('C2-A', [225, 87, 89])
        cmd.set_color('C2-B', [255, 157, 154])
        cmd.set_color('C3-A', [73, 152, 148])
        cmd.set_color('C3-B', [134, 188, 182])
        cmd.set_color('C4-A', [89, 161, 79])
        cmd.set_color('C4-B', [140, 209, 125])
        cmd.set_color('C5-A', [182, 153, 45])
        cmd.set_color('C5-B', [241, 206, 99])
        cmd.set_color('C6-A', [176, 122, 161])
        cmd.set_color('C6-B', [212, 166, 200])
        cmd.set_color('C7-A', [211, 114, 149])
        cmd.set_color('C7-B', [250, 191, 210])
        cmd.set_color('C8-A', [157, 118, 96])
        cmd.set_color('C8-B', [215, 181, 166])
        cmd.set_color('C9-A', [121, 112, 110])
        cmd.set_color('C9-B', [186, 176, 172])
    elif palette_name == 'plddt':
        for name, rgb, _ in _plldt_records:
            cmd.set_color(name, rgb)
    else:
        raise RuntimeError(f'Style {palette_name} not recognized')


cmd.extend(add_palette)


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
):
    """
    Load 3D arrow as object.

    Parameters
    ----------
    origin : list of three floats
        Vector of position of arrow tail.
    orient : list of three floats
        Vector of direction of arrow head. Does not need to be a unit vector.
    name : str
        Name of arrow object.
    color : str or list of three floats
        Registered color name or RGB values in range [0.0, 1.0] or [0, 255]
    tail_length : float
        Length of tail in A.
    tail_width : float
        Width of tail in A.
    head_length : float
        Length of head in A.
    head_width : float
        Width of head in A.
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

    cmd.load_cgo(arrow, name)


cmd.extend(load_cgo_arrow)


def load_nanobody_arrow(
    name,
    color,
    selection='*',
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
    color : str or list of three floats
        Registered color name or RGB values in range [0.0, 1.0] or [0, 255].
    selection : str
        Selection-expression or name-pattern corresponding to the nanobody atoms or object.
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
        raise RuntimeError(f'Mode {mode} not recognized. Must be in {accepted_modes}.')
    if arrow_kwargs is None:
        arrow_kwargs = {}

    # Gather CA atoms (with state assignments)
    atoms = []
    cmd.iterate_state(
        state,
        f'({selection}) and (name CA)',
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


cmd.extend(load_nanobody_arrow)


def load_nanobody_as_arrow(
    path,
    name=None,
    color=None,
    selection='*',
    arrow_kwargs=None,
    nanobody_kwargs=None,
):
    """
    Load a nanobody structure as an arrow.

    Parameters
    ----------
    path : str
        Path to structure.
    name : str
        Name of arrow object.
    color : str or list of three floats
        Registered color name or RGB values in range [0.0, 1.0] or [0, 255].
    selection : str
        Selection-expression or name-pattern corresponding to the nanobody atoms or object.
    arrow_kwargs : dict
        Additional arguments passed to load_nanobody_arrow call.
    nanobody_kwargs : dict
        Additional arguments unpacked in load_nanobody_arrow call.
    """
    if nanobody_kwargs is None:
        nanobody_kwargs = {}
    if arrow_kwargs is None:
        arrow_kwargs = {}

    # Some bookkeeping for internal clarity
    arrow_name = name
    del name

    cmd.load(path)
    model_name = cmd.get_object_list()[-1]

    # Set name and color
    if arrow_name is None:
        arrow_name = model_name
    if arrow_name == model_name:
        tmp_name = cmd.get_unused_name()
        cmd.set_name(model_name, tmp_name)
        model_name = tmp_name
    if color is None:
        indices = []
        cmd.iterate(
            f'({selection}) and %{model_name} and name CA', lambda x: indices.append(x.color)
        )
        index = indices[0]
        color = cmd.get_color_tuple(index)

    load_nanobody_arrow(
        arrow_name,
        color,
        f'({selection}) and %{model_name}',
        arrow_kwargs=arrow_kwargs,
        **nanobody_kwargs,
    )
    cmd.delete(f'%{model_name}')


cmd.extend(load_nanobody_as_arrow)
