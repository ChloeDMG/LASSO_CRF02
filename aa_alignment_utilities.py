#! /usr/bin/env python3

# Created on 10 apr. 2019
# Copyright (C) PTI CHU Purpan

__author__ = 'Nicolas JEANNE'
__copyright__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'jeanne.n@chu-toulouse.fr'

def prune_aln(alignment, list_idx):
    '''Remove the rows of an alignment using the index positions of a list.

    :param Bio.Align.MultipleSeqAlignment alignment: the alignment.
    :param list of int list_idx : list of index to remove
    :return: an alignment.
    :rtype: Bio.Align.MultipleSeqAlignment
    '''

    for index in sorted(list_idx, reverse=True):
        if index == 0:
            alignment = alignment[index+1:, :]
        elif index == len(alignment)-1:
            alignment = alignment[:index, :]
        else:
            tail_alignment = alignment[index+1:, :]
            alignment = alignment[:index, :]
            for i in range(len(tail_alignment)):
                alignment.append(tail_alignment[i, :])

    return alignment


def numbering_positions(aligned_reference):
    '''Numbering an alignment based on the aligned reference sequence.
    The positions are registred as posX. The gaps in the reference sequence
    which are insertions in the sequences of the alignement are registred as posX_1.
    IE:
        >reference_sequence
        CTR-PQ
        >foo
        CTRAPQ
        >bar
        C-KAPQ
    will produce:
        [pos1, pos2, pos3, pos3_1, pos4, pos5]

    :param str aligned_reference: the string of aligned reference sequence of an alignment.
    :return: the list of positions.
    :rtype: list of str
    '''

    positions_list = []
    position_no_gap = 0
    insertion = 1
    for aa in aligned_reference:
        if aa != '-':
            insertion = 1
            position_no_gap += 1
            positions_list.append('pos{}'.format(position_no_gap))
        else:
            positions_list.append('pos{}_{}'.format(position_no_gap, insertion))
            insertion += 1
    return positions_list


def numeric_recoding_from_seqrecord(sequence_record):
    '''Recode amino-acids in a Biopython sequence record with numerics.

    :param SeqRecord sequence_record: a Biopython amino acid sequence record object.
    :return: a list with the recoded sequence.
    :rtype: list of int.
    '''

    NUMERIC_AA = {'A': 1, 'R': 2, 'N': 3, 'D': 4, 'C': 5, 'E': 6, 'Q': 7, 'G': 8,
                  'H': 9, 'I': 10, 'L': 11, 'K': 12, 'M': 13, 'F': 14, 'P': 15,
                  'S': 16, 'T': 17, 'W': 18, 'Y': 19, 'V': 20, '-': 21, 'U': 22,
                  'O': 23}

    seq = sequence_record.seq
    recoded = []
    for aa in seq:
        recoded.append(NUMERIC_AA[aa])
    return recoded



def numeric_recoding_from_dataframe(dataframe, col_names_to_recode):
    '''Recode amino-acids in a dataframe with numerics.

    :param dataframe dataframe: the dataframe containing the amino-acids to recode.
    :param list col_names_to_recode: list of columns names to recode.
    :return: the recoded dataframe.
    :rtype: dataframe.
    '''

    NUMERIC_AA = {'A': 1, 'R': 2, 'N': 3, 'D': 4, 'C': 5, 'E': 6, 'Q': 7, 'G': 8,
                  'H': 9, 'I': 10, 'L': 11, 'K': 12, 'M': 13, 'F': 14, 'P': 15,
                  'S': 16, 'T': 17, 'W': 18, 'Y': 19, 'V': 20, '-': 21, 'U': 22,
                  'O': 23}

    # recode the selected columns
    tmp_df = dataframe[col_names_to_recode]
    for aa in NUMERIC_AA:
        tmp_df = tmp_df.replace(aa, NUMERIC_AA[aa])
    dataframe[col_names_to_recode] = tmp_df
    return dataframe
