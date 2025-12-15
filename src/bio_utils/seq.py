"""Functions for sequence data."""

from textwrap import dedent


def read_fasta(path):
    """
    Read FASTA file at path and return iterator of sequence records.

    See parse_fasta for details on parsing.

    Parameters
    ----------
    path : str
        Path to input FASTA file

    Returns
    -------
    iterator of tuples of (header, seq)
        The first element is the header line without the >, and the second element is the
        corresponding sequence.
    """
    with open(path) as file:
        yield from parse_fasta(file)


def write_fasta(path, records, width=100):
    """
    Write FASTA file at path using iterable of sequence records.

    Parameters
    ----------
    path : str
        Path to output FASTA file
    records : iterable of tuples of (header, seq)
        See read_fasta for details
    width : int
        Width of sequence lines
    """
    with open(path, 'w') as file:
        for header, seq in records:
            seqlines = '\n'.join(seq[i : i + width] for i in range(0, len(seq), width))
            file.write(f'>{header}\n{seqlines}\n')


def parse_fasta(file):
    """
    Read FASTA from stream and return iterator of sequence records.

    Parameters
    ----------
    file : file object
        File object of FASTA text stream

    Returns
    -------
    iterator of tuples of (header, seq)
    """
    line = file.readline()
    while line:
        if line.startswith('>'):
            header = line.rstrip()[1:]
            line = file.readline()

        seqlines = []
        while line and not line.startswith('>'):
            seqlines.append(line.rstrip())
            line = file.readline()
        seq = ''.join(seqlines)
        yield header, seq


def transcribe(dna_seq):
    """
    Transcribe a DNA sequence to RNA.

    Transcription is a simple mapping of T/t to U/u, respectively. No additional quality checks
    (e.g., for partially or already transcribed sequences) are performed.

    Parameters
    ----------
    dna_seq : str

    Returns
    -------
    rna_seq : str
    """
    table = {
        ord('T'): ord('U'),
        ord('t'): ord('u'),
    }
    rna_seq = dna_seq.translate(table)  # Confusing, but this is string translation
    return rna_seq


def translate(
    rna_seq,
    table=1,
    gap_syms=None,
    ambiguous_case='upper',
    strict_frame=True,
    strict_translation=True,
):
    """
    Translate an RNA sequence to protein.

    Translation is a simple mapping from codons to their corresponding amino acids. Stop codons,
    represented by an asterisk, are included in the output and "read-through," so an output sequence
    can contain multiple stop codons.

    Sequences of three gaps are converted to a single gap of the corresponding symbol. "Gap codons"
    with mixed symbols are converted to the first symbol in `gap_symbols`.

    The case of the output symbols matches their corresponding codons in the input. See
    `ambiguous_case` for mixed-case codons.

    Parameters
    ----------
    rna_seq : str
    table : int or dict
        If int, uses the corresponding NCBI translation table. If dict, provides a custom mapping of
        codons to amino acids. Custom mappings must provide codons in uppercase, as input
        sequences are converted to uppercase internally. Since non-alphabetic symbols are not
        uppercase or lowercase, the output case of codons with non-alphabetic symbols is controlled
        by `ambiguous_case`.
    gap_syms : str
        A string of symbols that are treated as gaps. Defaults to '-.'.
    ambiguous_case : {'upper', 'lower'}
        The case of output symbols for mixed-case codons.
    strict_frame : bool
        If True, raises a ValueError if the number of symbols is not a multiple of three. If False,
        incomplete codons are ignored.
    strict_translation : bool
        If True, raises a RuntimeError if a codon is not found in the translation table. If False,
        codons without a corresponding amino acid are converted to '?'.

    Returns
    -------
    protein_seq : str
    """
    if isinstance(table, int):
        table_string = table_strings[table]
        table = _parse_table(table_string)
    elif not isinstance(table, dict):
        raise ValueError('Parameter table is not int or dict')
    if gap_syms is None:
        gap_syms = '-.'
    elif not isinstance(gap_syms, str):
        raise ValueError('Parameter gap_syms is not str.')
    if ambiguous_case not in {'upper', 'lower'}:
        raise ValueError("Parameter ambiguous_case is not in {'upper', 'lower'}.")

    q, r = divmod(len(rna_seq), 3)
    if strict_frame and r != 0:
        raise ValueError('Number of symbols is not a multiple of three.')

    protein_seq = []
    for i in range(q):
        codon = rna_seq[3 * i : 3 * (i + 1)]
        all_gap = all(sym in gap_syms for sym in codon)
        if all_gap:
            if len(set(codon)) == 1:
                aa_sym = codon[0]
            else:
                aa_sym = gap_syms[0]
        else:
            all_upper = all(sym.isupper() for sym in codon)
            all_lower = all(sym.islower() for sym in codon)

            if codon.upper() in table:
                aa_sym = table[codon.upper()]
            elif strict_translation:
                raise RuntimeError(f'Codon {codon} not found in translation table.')
            else:
                aa_sym = '?'

            if all_upper:
                aa_sym = aa_sym.upper()
            elif all_lower:
                aa_sym = aa_sym.lower()
            elif ambiguous_case == 'upper':
                aa_sym = aa_sym.upper()
            elif ambiguous_case == 'lower':
                aa_sym = aa_sym.lower()
            else:
                raise RuntimeError("Parameter ambiguous case is not in {'upper', 'lower'}")
        protein_seq.append(aa_sym)
    protein_seq = ''.join(protein_seq)
    return protein_seq


def _parse_table(table_string):
    lines = table_string.split('\n')
    lines = [line.split('=')[1].strip() for line in lines]
    aa_syms = lines[0]
    b_lines = lines[2:5]
    b_lines = [line.replace('T', 'U') for line in b_lines]
    codons = [''.join(b_syms) for b_syms in zip(*b_lines)]
    table = {}
    for aa_sym, codon in zip(aa_syms, codons):
        table[codon] = aa_sym
    return table


table_strings = {
    1: dedent("""\
        AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
        Starts = ---M------**--*----M---------------M----------------------------
        Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
        Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
        Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"""),
    2: dedent("""\
        AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
        Starts = ----------**--------------------MMMM----------**---M------------
        Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
        Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
        Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"""),
    3: dedent("""\
        AAs  = FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
        Starts = ----------**----------------------MM---------------M------------
        Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
        Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
        Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"""),
}

# fmt: off
protein_alphabet = {
    'A', 'C', 'D', 'E',
    'F', 'G', 'H', 'I',
    'K', 'L', 'M', 'N',
    'P', 'Q', 'R', 'S',
    'T', 'V', 'W', 'Y',
    'X',
}
nucleic_alphabet = {
    'G', 'A', 'T', 'C',
    'U',
    'R', 'Y', 'M', 'K', 'S', 'W',
    'H', 'B', 'V', 'D',
    'N',
}
# fmt: on
