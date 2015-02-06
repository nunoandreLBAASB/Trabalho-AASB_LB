#!/usr/bin/env python3

from Bio import SeqIO
# from msvcrt import getch


def get_directory_filenames(directory):
    import os

    for _, _, files in os.walk(directory):
        for filename in files:
            yield filename


def get_uniprot_ids(query, organism):
    try:
        import urllib.request as urllib2
    except:
        import urllib2

    def ids(webpage, length=6):
        try:
            start = webpage.index('<tbody>')
        except ValueError:
            start = len(webpage)

        while start < len(webpage):
            try:
                start += webpage[start:].index('<tr id="') + len('<tr id="')
                # print(start, length, '|', webpage[start:length], '|', webpage[start:])
                # yield webpage[start:length]
                uid = ''
                for i in range(6):
                    uid += webpage[start + i]
                yield uid
            except ValueError:
                start = len(webpage)
        return

    uids = []
    query = query.replace(' ', '+')
    organism = organism.replace(' ', '+')
    url = 'http://www.uniprot.org/uniprot/?query="{0}"+AND+organism%3A"{1}"&sort=score'.format(query, organism)
    handle = urllib2.urlopen(url)
    webpage = str(handle.read())
    for uid in ids(webpage):
        uids.append(uid)
    return uids


def get_uniprot_id_data(uid, format):
    from Bio import SeqIO
    try:
        import urllib.request as urllib2
    except:
        import urllib2

    url = 'http://www.uniprot.org/uniprot/%s.%s' % (uid, format)
    page = urllib2.urlopen(url)
    return SeqIO.parse(page, 'uniprot-xml')


if __name__ == '__main__':
    directory = './ncbiSeqRecord'
    organism = 'Neisseria gonorrhoeae'

    for filename in get_directory_filenames(directory):
        with open(directory + '/' + filename, 'r') as input_handle:
            print()
            print('Fetching accession ID from file "%s"... ' % filename, end='', flush=True)
            accession = SeqIO.read(input_handle, 'gb').name
            print('Done')
            print(' ---\tAccession ID:', accession)
            print('Fetching Uniprot IDs... ', end='', flush=True)
            uids = get_uniprot_ids(accession, organism)
            print('Done')

            counter = 0
            filename = filename.split('.')[0]
            print(' ---\tUniprot IDs found:', ', '.join(uids))
            for uid in uids:
                filename = filename.split('_')[0] + ('_%d' % counter if counter else '')
                print('Fetching Uniprot ID "%s" data... ' % uid, end='', flush=True)
                record = get_uniprot_id_data(uid, 'xml')
                print('Done')

                print('Writing data to file...', end='', flush=True)
                # with open('uniprotSeqRecord/%s.xsd' % filename, 'w') as output_handle:
                #     SeqIO.write(record, output_handle, 'seqxml')
                with open('uniprotSeqRecord/%s.gb' % filename, 'w') as output_handle:
                    SeqIO.write(record, output_handle, 'gb')
                print('Done')

                counter += 1
