#!/usr/bin/env python3

from Bio.Blast import NCBIXML


def get_directory_filenames(directory):
    import os

    for _, _, files in os.walk(directory):
        for filename in files:
            yield filename


if __name__ == '__main__':
    input_dir = '../results/blast/swissprot'
    output_file = '../results/functions.txt'

    with open(output_file, 'w') as f:
        pass  # create/clear file

    for filename in get_directory_filenames(input_dir):
        print(filename)
        with open('%s/%s' % (input_dir, filename), 'r') as blast_file:
            blast_records = NCBIXML.parse(blast_file)
            for rec in blast_records:
                saved_hsp = None
                for alignment in rec.alignments:
                    for hsp in alignment.hsps:
                        try:
                            if saved_hsp.expect > hsp.expect:
                                saved_alignment = alignment
                                saved_hsp = hsp
                        except:
                            saved_alignment = alignment
                            saved_hsp = hsp
                # print('\n****Alignment****')
                # print('sequence:', saved_alignment.title)
                # print('length:', saved_alignment.length)
                # print('def:', saved_alignment.hit_def)
                # print('id:', saved_alignment.hit_id)
                # print('strand:', saved_hsp.strand)
                # print('e value:', saved_hsp.expect)
                # print('identity:', saved_hsp.score)
                # print(saved_hsp.bits)
                # print(saved_hsp.query[0:75] + '...')
                # print(saved_hsp.match[0:75] + '...')
                # print(saved_hsp.sbjct[0:75] + '...')
            if saved_hsp:
                ngo = filename.split('.')[0]
                function = saved_alignment.hit_def
                start = function.index('=') + 1
                try:
                    end = function.index(';')
                except ValueError:
                    try:
                        end = function.index('[')
                    except ValueError:
                        end = None
                with open(output_file, 'a') as f:
                    f.write('%s\t%s\t%s\tSimilar to: %s\n' % (ngo, saved_hsp.expect, saved_hsp.score, function))
