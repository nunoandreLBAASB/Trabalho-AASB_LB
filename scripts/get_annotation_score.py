#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
# from msvcrt import getch


def get_directory_filenames(directory):
    import os

    for _, _, files in os.walk(directory):
        for filename in files:
            yield filename


def get_uniprot_id_data(uid, data_format=''):
    from Bio import SeqIO
    try:
        import urllib.request as urllib2
    except:
        import urllib2

    if data_format:
        data_format = '.' + data_format
    url = 'http://www.uniprot.org/uniprot/%s%s' % (uid, data_format)
    page = urllib2.urlopen(url)
    return str(page.read())


def get_annotation_status(uniprot_webpage):
    status = []

    try:
        start = uniprot_webpage.index('Unreviewed')
        reviewed = '0'
    except ValueError:
        try:
            start = uniprot_webpage.index('Reviewed')
            reviewed = '1'
        except ValueError:
            reviewed = 'ND'
    status.append(reviewed)

    try:
        start += uniprot_webpage[start:].index('toolTipContent')
        start += uniprot_webpage[start:].index('Annotation score: ') + len('Annotation score: ')
        score = uniprot_webpage[start]
    except ValueError:
        score = 'ND'
    status.append(score)

    try:
        start += uniprot_webpage[start:].index('sep')
        start += uniprot_webpage[start:].index('</span>') + len('</span>')
        end = start + uniprot_webpage[start:].index('<')
        description = uniprot_webpage[start:end]
    except ValueError:
        description = 'ND'
    status.append(description)

    return status


if __name__ == '__main__':
    input_dir = '../data/uniprotSequence'
    output_file = '../results/score.txt'
    extension = 'gb'
    organism = 'Neisseria gonorrhoeae'
    info = ''

    for filename in get_directory_filenames(input_dir):
        line = ''
        if filename.split('.')[-1] == extension:
            line = filename.split('.')[0] + '\t'

            with open(input_dir + '/' + filename, 'r') as input_handle:
                accession = SeqIO.read(input_handle, extension).id
            uniprot_record = get_uniprot_id_data(accession)

            status = get_annotation_status(uniprot_record)
            line += '\t'.join(status)
            print(line)
        info += line + '\n'

    with open(output_file, 'w') as f:
        f.write(info)
