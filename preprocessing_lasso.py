#! /usr/bin/env python3

import argparse
import os
import csv
import re
import shutil
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import aa_alignment_utilities as util_aln
import aa_sequences_utilities as util_seq


__author__ = 'Nicolas JEANNE'
__copyright__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'jeanne.n@chu-toulouse.fr'

'''
Created on 11 jun. 2019

# Copyright (C) 2019 PTI CHU Purpan
'''

if __name__ == '__main__':
    # Argument parsing
    descr = '''
    {} v. {}

    Created by {}.
    Contact: {}
    {}

    Distributed on an 'AS IS' basis without warranties or conditions of any kind, either express or implied.

    Pre-processing of the fasta AA sequences for the MatLab LASSO step.

    Align the AA fasta sequences of HIV-1 enveloppe provided as input with the reference sequences.
    Add the charge and N-glycosylation sites infos for the provided sequence and
    convert the alignment to a CSV file file which amino acids are recoded numerically.

    The input fasta file must contains HIV-1 envelope sequences to predict.

    The output directory will contain sub-directories with the ID <ID_SEQ> of each sequence to predict. The subdirectories contains
    a CSV file with the LASSO training sequences <ID_SEQ>_to_train.csv and a CSV file of the sequence to predict <ID_SEQ>_to_predict.csv.
    '''.format(os.path.splitext(os.path.basename(__file__))[0], __version__, __author__, __email__, __copyright__)
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--out', required=True, help='path to the output directory.')
    parser.add_argument('input', help='path to the fasta file. The file contains unaligned amino acid sequences of HIV-1 env of which phenotropism have to be predicted.')
    args = parser.parse_args()

    # create result directory
    out_dir = os.path.abspath(args.out)
    tmp = os.path.join(out_dir, 'tmp')
    if not os.path.exists(tmp):
        os.makedirs(tmp, exist_ok=True)

    # load the reference AA sequences
    ref_seq_path = os.path.join(os.path.dirname(__file__), 'ref', 'CRF02_ref_seq.fa')

    references = list(SeqIO.parse(ref_seq_path, 'fasta'))

    # load the input AA sequence
    input_records = SeqIO.parse(args.input, 'fasta')
    for seq_record in input_records:
        seq_to_predict_id = seq_record.id
        print('Processing: {}\n\tAlignment with Muscle, please wait..'.format(seq_to_predict_id))
        # create the output directory for the sequence to predict
        res_dir_path = os.path.join(out_dir, seq_to_predict_id)
        os.makedirs(res_dir_path, exist_ok=True)

        # copy the references sequences used in LASSO
        fasta_to_align = references.copy()
        # add the sequence to predict
        fasta_to_align.append(seq_record)

        # align with Muscle the AA sequences and the reference sequence
        fasta_to_align_path = os.path.join(tmp, '{}_to_align.fa'.format(seq_to_predict_id))
        SeqIO.write(fasta_to_align, fasta_to_align_path, 'fasta')
        aln_path = os.path.join(tmp, '{}_aln.fa'.format(seq_to_predict_id))
        muscle_cmd = MuscleCommandline(input=fasta_to_align_path, out=aln_path)
        stdout, stderr = muscle_cmd()

        # cut the alignment on the reference length
        alignment = AlignIO.read(aln_path, 'fasta')
        cut_start = None
        cut_end = None
        for i, seqrecord in enumerate(alignment):
            if seqrecord.id == 'reference':
                aln_ref_seq = seqrecord.seq
                cut_start = 0
                cut_end = len(aln_ref_seq) - 1
                # look for the first position which is not a gap in the alignment
                while aln_ref_seq[cut_start] == '-':
                    cut_start += 1
                # look for the first position which is not a gap in the alignment starting from the end
                while aln_ref_seq[cut_end] == '-':
                    cut_end -= 1
                break
        # cut on the positions which are not only gaps at the begining and the end of the alignment
        AlignIO.write(alignment[:, cut_start:cut_end+1], aln_path, 'fasta')

        # extract info from the alignment
        aln_pos_numbering = None
        data = {}
        pattern_id = re.compile('(R5|X4|R5X4)_([A-Za-z0-9-]+)-(\\d+)$')
        for aln in SeqIO.parse(aln_path, 'fasta'):
            if aln.id == 'reference':
                # numbering positions of the alignment based on the reference sequence
                aln_pos_numbering = util_aln.numbering_positions(aln)
            else:
                if aln.id != seq_to_predict_id:
                    searched_result = pattern_id.search(aln.id)
                    if searched_result:
                        sample = searched_result.group(2)
                        clone = searched_result.group(3)
                        tropism = searched_result.group(1)
                else:
                    sample = seq_to_predict_id
                    clone = 'N/A'
                    tropism = None
                # recode the amino acids with numerics
                recoded_aa = util_aln.numeric_recoding_from_seqrecord(aln)
                # compute the charge of the sequence
                charge = util_seq.net_charge(aln)
                # compute the number of N-glycosylation sites
                nglyco = util_seq.count_N_glycosylation_site(aln)[0]
                # record the data
                row = [sample, clone]
                for item in recoded_aa:
                    row.append(item)
                row.append(charge)
                row.append(nglyco)
                data[aln.id] = {'row': row,
                                'tropism': tropism}


        # create the CSV for the sequence to predict
        to_predict_path = os.path.join(res_dir_path, '{}_to_predict.csv'.format(seq_to_predict_id))
        print('\tWriting the CSV file for prediction: {}'.format(to_predict_path))
        header = ['sample', 'clone'] + aln_pos_numbering + ['net_charge', 'Nglyco']
        data_to_predict = data.pop(seq_to_predict_id)
        with open(to_predict_path, 'w') as predict_csv:
            writer = csv.writer(predict_csv, delimiter=';')
            writer.writerow(header)
            writer.writerow(data_to_predict)

        # create the CSV for the sequences used for the LASSO training
        to_train_path = os.path.join(res_dir_path, '{}_to_train.csv'.format(seq_to_predict_id))
        print('\tWriting the CSV file for training: {}'.format(to_predict_path))
        header.append('predictor')
        with open(to_train_path, 'w') as train_csv:
            writer = csv.writer(train_csv, delimiter=';')
            writer.writerow(header)
            for seq_id in data:
                row = data[seq_id]['row']
                row.append(data[seq_id]['tropism'])
                writer.writerow(row)

    shutil.rmtree(tmp)
