#!/usr/bin/env python3

from Bio import SeqIO
# from msvcrt import getch


def get_directory_filenames(directory):
    import os

    for _, _, files in os.walk(directory):
        for filename in files:
            yield filename


def get_uniprot_id_data(uid, format):
    from Bio import SeqIO
    try:
        import urllib.request as urllib2
    except:
        import urllib2

    url = 'http://www.uniprot.org/uniprot/%s.%s' % (uid, format)
    page = urllib2.urlopen(url)
    return str(page.read())


def is_reviewed(rdf):
    boolean = {'true': True, 'false': False}
    tag = 'reviewed'
    s, e = '>', '<'
    try:
        start = rdf.index(tag)
        start = start + rdf[start:].index(s) + 1
        end = start + rdf[start:].index(e)
        value = rdf[start:end].lower()
        value = ''
        for i in range(end - start):
            value += rdf[start + i]
        return boolean[value.lower()]
    except ValueError:
        return False


if __name__ == '__main__':
    directory = './uniprotSeqRecord'
    organism = 'Neisseria gonorrhoeae'
    notreviewed = ''

    for filename in get_directory_filenames(directory):
        if filename.split('.')[-1] == 'gb':
            with open(directory + '/' + filename, 'r') as input_handle:
                print(filename.split('.')[0], '->', end=' ', flush=True)
                accession = SeqIO.read(input_handle, 'gb').id
                uniprot_record = get_uniprot_id_data(accession, 'rdf')
                if not is_reviewed(uniprot_record):
                    print('Not reviewed.')
                    notreviewed += filename.split('.')[0] + '\n'
                else:
                    print('Reviewed!')

    with open('notreviewed.txt', 'w') as f:
        f.write(notreviewed)
