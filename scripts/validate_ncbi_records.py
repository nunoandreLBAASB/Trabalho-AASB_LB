#!/usr/bin/env python3


# verify if information in the feature is the same as the one present in the line
def verify(line, feature, ltstart, ltend):
    check = False
    start, end, strand = feature.location.start + 1, feature.location.end, feature.location.strand
    if start == int(line[2]) and end == int(line[3]) and strand == int(line[4] + '1'):
        check = True
    return check

from Bio import Entrez, SeqIO
Entrez.email = "andremsantiago_@hotmail.com"
gi = '59717368'
ltstart, ltend = 'NGO0971', 'NGO1212'


def import_seqrecord(gi):
    # fetches the full genome and converts to a SeqRecord
    request = Entrez.efetch(db="nucleotide", rettype='gb', retmode='text', id=gi)
    record = SeqIO.read(request, 'gb')
    request.close()
    return record


def read_seqrecord(filename, format):
    request = open(filename, 'r')
    record = SeqIO.read(request, format)
    request.close()
    return record


record = read_seqrecord('./fullSeq/sequence_full.gb', 'gb')

# open log file to record validation
with open('log.log', 'w') as f:

    # open comparison table
    with open('ProteinTable864_169534.tsv') as table:
        for line in table:
            if line[0] != '#':
                line = line[:-1].split('\t')

                # if current locus_tag is between target locus_tag's
                if line[7] >= ltstart and line[7] <= ltend:
                    ngos = []

                    # fetch features with identical locus_tag
                    for feature in record.features:
                        if feature.type == 'gene' or feature.type == 'CDS':
                            try:
                                if feature.qualifiers['locus_tag'][0] == line[7]:
                                    ngos.append(feature)
                            except:
                                pass

                    # compare if information in each feature is the same as the one present in the comparison table and write to file
                    for ngo in ngos:
                        if verify(line, ngo, ltstart, ltend):
                            f.write('Check! ' + line[7] + ' ' + str(ngo.qualifiers['locus_tag'][0]) + ' ' + str(ngo.type) + '\n')
                        else:
                            f.write('Not check...\n' + str(ngo) + str(line) + '\n')
