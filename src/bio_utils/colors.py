"""Colors and functions for manipulating them."""

from string import hexdigits

# COLOR SCHEMES
color_schemes = {}

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
color_schemes['tableau_10'] = tableau_10

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
color_schemes['tableau_20'] = tableau_20

tableau_blacks = [
    '#000000',
    '#1B1B1B',
    '#333333',
    '#555555',
    '#666666',
    '#898989',
]
color_schemes['tableau_blacks'] = tableau_blacks

tableau_whites = [
    '#FFFFFF',
    '#F5F5F5',
    '#E6E6E6',
    '#D4D4D4',
    '#C0C0C0',
    '#B4B4B4',
]
color_schemes['tableau_whites'] = tableau_whites

tableau_oranges = [
    '#F28E2B',
    '#F9A655',
    '#FFBE7D',
]
color_schemes['tableau_oranges'] = tableau_oranges

tableau_yellows = [
    '#B6992D',
    '#D3B348',
    '#F1CE63',
]
color_schemes['tableau_yellows'] = tableau_yellows

tableau_greens = [
    '#59A14F',
    '#72B966',
    '#8CD17D',
]
color_schemes['tableau_greens'] = tableau_greens

tableau_cyans = [
    '#499894',
    '#69AAA5',
    '#86BCB6',
]
color_schemes['tableau_cyans'] = tableau_cyans

tableau_blues = [
    '#4E79A7',
    '#75A1C7',
    '#A0CBE8',
]
color_schemes['tableau_blues'] = tableau_blues

tableau_blue_greys = [
    '#606B76',
    '#818B91',
    '#A7ACAD',
]
color_schemes['tableau_blue_greys'] = tableau_blue_greys

tableau_purples = [
    '#B07AA1',
    '#C290B4',
    '#D4A6C8',
]
color_schemes['tableau_purples'] = tableau_purples

tableau_pinks = [
    '#D37295',
    '#E799B3',
    '#FABFD2',
]
color_schemes['tableau_pinks'] = tableau_pinks

tableau_reds = [
    '#E15759',
    '#F17B79',
    '#FF9D9A',
]
color_schemes['tableau_reds'] = tableau_reds

tableau_browns = [
    '#9D7660',
    '#BA9582',
    '#D7B5A6',
]
color_schemes['tableau_whites'] = tableau_whites

tableau_greys = [
    '#79706E',
    '#998F8C',
    '#BAB0AC',
]
color_schemes['tableau_browns'] = tableau_browns

tableau_highlights = [
    '#76B7B2',
    '#FF9DA7',
    '#EDC948',
]
color_schemes['tableau_highlights'] = tableau_highlights

# COLOR MAPS
color_maps = {}

protein_color_map = {
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
nucleic_color_map = {
    'A': '#6DD7A1',
    'T': '#EE8485',
    'G': '#E2DEDD',
    'C': '#7FADEA',
    'N': '#93908F',
}
