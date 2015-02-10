#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from sys import stdout

def get_directory_filenames(directory):
    import os

    for _, _, files in os.walk(directory):
        for filename in files:
            yield filename


def set_directory(directory):
    import os

    if not os.path.exists(directory):
        os.makedirs(directory)
        return True
    return False


if __name__ == '__main__':
    input_directory = '../data/ncbiSequence'
    output_directory = '../results/blast/ncbi'
    score = '../results/score.txt'
    extension = 'gb'
    database = 'nr' # 'swissprot'  # 'nr'
    program = 'blastp'
    replace = False
    notreviewed = []
    # E_VALUE_THRESH = 0.04

    set_directory('%s/%s' % (output_directory, database))

    with open(score, 'r') as f:
        scores = f.readlines()
    for line in scores:
        if len(line) > 1:
            line = line.split('\t')
            ngo = line[0]
            reviewed = bool(int(line[1]))
            if not reviewed:
                notreviewed.append(ngo)

    existent = []
    if not replace:
        for filename in get_directory_filenames(output_directory):
            filename = filename.split('.')
            filename[1] = 'gb'
            filename = '.'.join(filename)
            existent.append(filename)

    for filename in get_directory_filenames(input_directory):
        # if filename not in existent and filename.split('.')[0][-1] in notreviewed and filename.split('.')[1] == extension:
        if filename.split('.')[0][-1] == 'a':
            print()
            print('Fetching sequence record from "%s"... ' % filename)
            stdout.flush()
            with open(input_directory + '/' + filename, 'r') as input_handle:
                record = SeqIO.read(input_handle, format='gb')
            print('Done')
            filename = filename.split('.')[0]
            print('Perfoming %s against "%s" database... ' % (program, database))
            stdout.flush()
            result_handle = NCBIWWW.qblast(program, database, record.format('fasta'))
            print('Done')
            with open('%s/%s/%s.xml' % (output_directory, database, filename), 'w') as f:
                f.write(result_handle.read())
