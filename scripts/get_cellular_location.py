#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
# from msvcrt import getch


def get_directory_filenames(directory):
    import os

    for _, _, files in os.walk(directory):
        for filename in files:
            yield filename


def get_uniprot_id_data(uid, format=''):
    from Bio import SeqIO
    try:
        import urllib.request as urllib2
    except:
        import urllib2

    if format:
        format = '.' + format
    url = 'http://www.uniprot.org/uniprot/%s%s' % (uid, format)
    page = urllib2.urlopen(url)
    return str(page.read())


def get_cellular_locations(page):
    info = []
    start = 0

    while start < len(page):
        try:
            start += page[start:].index('/locations/SL-') + len('/locations/SL-')
            location_id = page[start:start + 4]

            start += page[start:].index('>') + 1
            end = start + page[start:].index('<')
            location_name = page[start:end]

            info.append([location_id, location_name])
        except ValueError:
            start = len(page)

    return info


if __name__ == '__main__':
    input_dir = '../data/uniprotSequence'
    output_file = '../results/location.txt'
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

            location = get_cellular_locations(uniprot_record)
            for i in range(len(location)):
                location[i] = '\t'.join(location[i])
            line += '\n\t'.join(location)
            print(line)
        info += line + '\n'

    with open(output_file, 'w') as f:
        f.write(info)
