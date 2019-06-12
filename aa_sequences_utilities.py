#! /usr/bin/env python3

# Created on 10 apr. 2019
# Copyright (C) PTI CHU Purpan

import sys
import traceback
import re

__author__ = 'Nicolas JEANNE'
__copyright__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'jeanne.n@chu-toulouse.fr'


def check_if_aa_seq(assumed_aa_seq, id_seq):
    '''Check if the SeqRecord is a valid amino-acids sequence.

    :param SeqRecord assumed_aa_record: the BioPython SeqRecord object.
    '''
    valid_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
                'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'U', 'O']
    if False in [aa in valid_aa for aa in assumed_aa_seq]:
        raise ValueError('{}: the amino-acids sequence contains at least an unknown character not included in:\n{}'.format(id_seq, ', '.join(valid_aa)))


def count_overlapping(pattern_string, string_to_search):
    ''''
    Find occurrence count of overlapping pattern and get their count and index positions.
    Start from left and start searching for the pattern, when found increment the counter
    and keep on search from next index position.

    :param str pattern_string: the regular expression string.
    :param str string_to_search: the string in which the pattern is searched.
    :return: a tuple of the number of occurences and the list of their start positions.
    :rtype: tuple(int, list)
    '''
    counter = 0
    start = 0
    indexes = []
    pattern = re.compile(pattern_string)
    while True:
        match = pattern.search(string_to_search, start)
        if match is None:
            return counter, indexes
        indexes.append(match.start())
        counter += 1
        start = 1 + match.start()

def net_charge(aa_record, use_H=False):
    '''Compute the net charge of an amino-acids sequence (BioPython SeqRecord) using:
        R, K and optionaly H residues: +1
        D, E residues: -1

    :param Bio.SeqRecord aa_record: the biopython amino-acids sequence record.
    :param boolean use_H: use H residues as +1 for the net charge computation, default is False.
    :return: the net charge of the amin acid sequence.
    :rtype: int
    '''

    try:
        if not isinstance(use_H, bool):
            raise TypeError('{}The type of use_H argument is not a boolean.'.format(traceback.format_exc()))
        # remove the gaps from the sequence and check if it is a valid AA sequence
        seq = aa_record.seq.ungap('-').upper()
        seq_id = aa_record.id
        check_if_aa_seq(seq, seq_id)
    except ValueError:
        print(traceback.format_exc())
        sys.exit(1)
    except TypeError:
        print(traceback.format_exc())
        sys.exit(1)

    # compute the net charge
    computed_charge = 0
    if use_H:
        positive_aa = ['R', 'K', 'H']
    else:
        positive_aa = ['R', 'K']
    negative_aa = ['D', 'E']
    for aa in seq:
        if aa in positive_aa:
            computed_charge += 1
        elif aa in negative_aa:
            computed_charge -= 1

    return computed_charge


def net_charge_pk(aa_record):
    '''Compute the net charge of an amino-acids sequence (BioPython SeqRecord) using
    the dissociation constants of the amino-acids and the physiological pH (7.4).

    :param SeqRecord aa_record: the biopython amino-acids sequence record.
    :return: the net charge of the amin acid sequence.
    :rtype: float
    '''

    # Methode SeqinR
    # REF: https://www.rdocumentation.org/packages/seqinr/versions/3.1-3/topics/pK
    print('''[WARNING] net_charge_pk(), charge Net calculée à partir de la méthode SeqinR (https://www.rdocumentation.org/packages/seqinr/versions/3.1-3/topics/pK)
          Différente de la méthode MOORE net_charge_pk2(), reste à déterminer laquelle est la plus adaptée.''')

    try:
        # remove the gaps from the sequence and check if it is a valid AA sequence
        seq = aa_record.seq.ungap('-').upper()
        seq_id = aa_record.id
        check_if_aa_seq(seq, seq_id)
    except ValueError:
        print(traceback.format_exc())
        sys.exit(1)
    except TypeError:
        print(traceback.format_exc())
        sys.exit(1)

    # compute the net charge
    pH = 7.4
    positive_aa = ['R', 'K', 'H']
    negative_aa = ['D', 'E', 'C', 'Y']

    pK = {'NH3': {'A': 9.69, 'R': 9.04, 'N': 8.8, 'D': 9.6, 'C': 10.28, 'E': 9.67,
                  'Q': 9.13, 'G': 9.6, 'H': 9.17, 'O': 9.65, 'I': 9.6, 'L': 9.6,
                  'K': 8.95, 'M': 9.21, 'F': 9.13, 'P': 10.6, 'S': 9.15, 'T': 9.1,
                  'W': 9.39, 'Y': 9.11, 'V': 9.62, 'U': 0.0},
          'COOH': {'A': 2.34, 'R': 2.17, 'N': 2.02, 'D': 1.88, 'C': 1.96, 'E': 2.19,
                   'Q': 2.17, 'G': 2.34, 'H': 1.82, 'O': 1.82, 'I': 2.36, 'L': 2.36,
                   'K': 2.18, 'M': 2.28, 'F': 1.83, 'P': 1.99, 'S': 2.21, 'T': 2.09,
                   'W': 2.83, 'Y': 2.2, 'V': 2.32, 'U': 0.0},
          'lateral': {'R' :  12.48, 'H' :  6, 'K' :  10.53, 'D' :  3.65, 'E' :  4.25, 'C': 8.18, 'Y': 10.07}}

    # get the N-Terminal AA pK, if it is a X (unknown AA from ambiguous codon), pK = 0 otherwhise compute it.
    if aa_record.seq[0] != 'X':
        pK_NH3 = pK['NH3'][aa_record.seq[0]]
    else:
        pK_NH3 = 0.0

    # compute the positive charge
    positive_charge = 10.0**(-pH) / (10.0**(-pK_NH3) + 10.0**(-pH))
    # for each positive AA get the number of ocurrences and its pK charge, then
    # compute the charge for this AA and add it to the positive charge
    for aa in positive_aa:
        nb_aa = len([aa for aa in seq[1:-1] if aa in positive_aa])
        pK_lateral = pK['lateral'][aa]
        positive_charge += nb_aa * 10.0**(-pH) / (10.0**(-pK_lateral) + 10.0**(-pH))

    # get the C-Terminal AA), if it is a X (unknown AA from ambiguous codon), pK_COOH = 0 else compute it.
    if aa_record.seq[-1] != 'X':
        pK_COOH = pK['COOH'][aa_record.seq[-1]]
    else:
        pK_COOH = 0.0

    # compute the negative charge
    negative_charge = 10.0**(-pK_COOH) / (10.0**(-pK_COOH) + 10.0**(-pH))
    # for each negative AA get the number of ocurrences and its pK charge, then
    # compute the charge for this AA and add it to the negative charge
    for aa in negative_aa:
        nb_aa = len([aa for aa in seq[1:-1] if aa in negative_aa])
        pK_lateral = pK['lateral'][aa]
        negative_charge += nb_aa * 10.0**(-pK_lateral) / (10.0**(-pK_lateral) + 10.0**(-pH))

    # compute the final charge
    computed_charge = round(positive_charge - negative_charge, 2)
    return computed_charge

def net_charge_pk2(aa_record):
    '''Compute the net charge of an amino-acids sequence (BioPython SeqRecord) using
    the dissociation constants of the amino-acids and the physiological pH (7.4).

    :param SeqRecord aa_record: the biopython amino-acids sequence record.
    :return: the net charge of the amin acid sequence.
    :rtype: float
    '''

    # Methode MOORE
    # REF: https://onlinelibrary.wiley.com/doi/epdf/10.1016/0307-4412%2885%2990114-1
    print('''[WARNING] net_charge_pk2(), charge Net calculée à partir de la méthode MOORE (https://onlinelibrary.wiley.com/doi/epdf/10.1016/0307-4412%2885%2990114-1)
          Différente de la méthode SeqinR net_charge_pk(), reste à déterminer laquelle est la plus adaptée.''')


    try:
        # remove the gaps from the sequence and check if it is a valid AA sequence
        seq = aa_record.seq.ungap('-').upper()
        seq_id = aa_record.id
        check_if_aa_seq(seq, seq_id)
    except ValueError:
        print(traceback.format_exc())
        sys.exit(1)
    except TypeError:
        print(traceback.format_exc())
        sys.exit(1)

    # compute the net charge
    pH = 7.4
    positive_aa = ['R', 'K', 'H']
    negative_aa = ['D', 'E']

    pK = {'NH3': {'A': 9.69, 'R': 9.04, 'N': 8.8, 'D': 9.6, 'C': 10.28, 'E': 9.67,
                  'Q': 9.13, 'G': 9.6, 'H': 9.17, 'O': 9.65, 'I': 9.6, 'L': 9.6,
                  'K': 8.95, 'M': 9.21, 'F': 9.13, 'P': 10.6, 'S': 9.15, 'T': 9.1,
                  'W': 9.39, 'Y': 9.11, 'V': 9.62, 'U': 0.0},
          'COOH': {'A': 2.34, 'R': 2.17, 'N': 2.02, 'D': 1.88, 'C': 1.96, 'E': 2.19,
                   'Q': 2.17, 'G': 2.34, 'H': 1.82, 'O': 1.82, 'I': 2.36, 'L': 2.36,
                   'K': 2.18, 'M': 2.28, 'F': 1.83, 'P': 1.99, 'S': 2.21, 'T': 2.09,
                   'W': 2.83, 'Y': 2.2, 'V': 2.32, 'U': 0.0},
          'lateral': {'R' :  12.48, 'H' :  6, 'K' :  10.53, 'D' :  3.65, 'E' :  4.25}}

    # get the N-Terminal AA pK, if it is a X (unknown AA from ambiguous codon), pK = 0 otherwhise compute it.
    if aa_record.seq[0] != 'X':
        pK_NH3 = pK['NH3'][aa_record.seq[0]]
    else:
        pK_NH3 = 0.0

    # compute the positive charge
    positive_charge = 1.0 / (1.0 + 10.0**(pH - pK_NH3))
    # for each positive AA get the number of ocurrences and its pK charge, then
    # compute the charge for this AA and add it to the positive charge
    for aa in positive_aa:
        nb_aa = len([aa for aa in seq[1:-1] if aa in positive_aa])
        pK_lateral = pK['lateral'][aa]
        positive_charge += nb_aa / (1.0 + 10.0**(pH - pK_lateral))

    # get the C-Terminal AA), if it is a X (unknown AA from ambiguous codon), pK_COOH = 0 else compute it.
    if aa_record.seq[-1] != 'X':
        pK_COOH = pK['COOH'][aa_record.seq[-1]]
    else:
        pK_COOH = 0.0

    # compute the negative charge
    negative_charge = -1.0 / (1.0 + 10.0**(-pH + pK_COOH))
    # for each negative AA get the number of ocurrences and its pK charge, then
    # compute the charge for this AA and add it to the negative charge
    for aa in negative_aa:
        nb_aa = len([aa for aa in seq[1:-1] if aa in negative_aa])
        pK_lateral = pK['lateral'][aa]
        negative_charge += -nb_aa / (1.0 + 10.0**(-pH + pK_lateral))

    # compute the final charge
    computed_charge = round(positive_charge - negative_charge, 2)
    return computed_charge

def count_N_glycosylation_site(aa_record):
    '''Count the number of N glycosilation sites in the amino-acids sequence using
    the pattern N[^P][TS][^P] defined from https://www.hiv.lanl.gov/content/sequence/GLYCOSITE/glycosite.html

    :param SeqRecord aa_record: the biopython amino-acids sequence record.
    :return: a tuple of the number N-glycosylation sites and the list of their start positions
    :rtype: tuple(int, list)
    '''
    try:
        # remove the gaps from the sequence and check if it is a valid AA sequence
        seq = aa_record.seq.ungap('-').upper()
        seq_id = aa_record.id
        check_if_aa_seq(seq, seq_id)
    except ValueError:
        print(traceback.format_exc())
        sys.exit(1)
    except TypeError:
        print(traceback.format_exc())
        sys.exit(1)

    pattern_nglyco = 'N[^P][TS][^P]'
    return count_overlapping(pattern_nglyco, str(seq))
