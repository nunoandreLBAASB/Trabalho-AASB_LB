#!/usr/bin/env python3

from Bio import Entrez, SeqIO
Entrez.email = 'andremsantiago_@hotmail.com'


def import_seqrecord(gi):
    # fetches the full genome and converts to a SeqRecord
    request = Entrez.efetch(db='gene', rettype='docsum', retmode='text', id=gi)
    record = request.read()
    request.close()
    return record


if __name__ == '__main__':
    ltstart, ltend = 'NGO0971', 'NGO1212'
    input_file = '../data/ncbiGenome/sequence_full.gb'
    output_dir = '../data/ncbiGene'
    with open(input_file, 'r') as request:
        record = SeqIO.read(request, 'gb')
        for feature in record.features:
            if feature.type == 'gene':
                try:
                    print(feature.qualifiers['locus_tag'][0])
                    if ltend >= feature.qualifiers['locus_tag'][0] >= ltstart:
                        geneid = feature.qualifiers['db_xref']
                        for gi in geneid:
                            if gi.split(':')[0] == 'GeneID':
                                geneid = gi.split(':')[1]
                                break
                        with open('%s/%s.txt' % (output_dir, feature.qualifiers['locus_tag'][0]), 'w') as f:
                            f.write(import_seqrecord(geneid))
                except:
                    pass
