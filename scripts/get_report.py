#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Blast import NCBIXML


def get_refids(filename):
    refids = {}
    with open(filename, 'r') as f:
        for line in f:
            if line[:11] == 'Search term':
                start = line.index('AND ') + len('AND ')
                end = start + line[start:].index('[')
                product = line[start:end]
                refids[product] = []
            elif line[:3] == 'IDs':
                line = line.split(' ')[1:]
                refids[product].extend(line)
    return refids


def get_location(filename):
    location = {}
    with open(filename, 'r') as f:
        for line in f:
            if len(line) > 1:
                line = line.rstrip().split('\t')
                if line[0]:
                    locus_tag = line[0]
                if line[1:]:
                    if locus_tag not in location:
                        location[locus_tag] = []
                    location[locus_tag].append([line[1], line[2]])
    return location


def get_anotation_score(filename):
    score = {}
    with open(filename, 'r') as f:
        for line in f:
            if len(line) > 1:
                line = line.rstrip().split('\t')
                score[line[0]] = [line[1], line[2], line[3]]
    return score


def get_blast_results(filename):
    blast = {}
    with open(filename, 'r') as f:
        blast_records = NCBIXML.parse(f)
        for record in blast_records:
            for alignment in record.alignments:
                blast[alignment.hit_id] = []
                blast[alignment.hit_id].append(alignment.hit_def)
                blast[alignment.hit_id].append(alignment.length)
                for hsp in alignment.hsps:
                    blast[alignment.hit_id].append(hsp.expect)
                    blast[alignment.hit_id].append(hsp.score)
                    blast[alignment.hit_id].append(hsp.query)
                    blast[alignment.hit_id].append(hsp.match)
                    blast[alignment.hit_id].append(hsp.sbjct)
    return blast


if __name__ == '__main__':
    output_file = '../results/report.txt'
    with open(output_file, 'w') as f:
        pass
    with open('../results/ngos.txt', 'r') as f:
        ngos = f.read()
    ngos = ngos.rstrip().split('\n')
    annotation_scores = get_anotation_score('../results/score.txt')
    locations = get_location('../results/location.txt')
    ref_ids = get_refids('../results/ref_ids.txt')

    for ngo in ngos:
        print(ngo)
        ncbi, uniprot = None, None
        try:
            blasts = get_blast_results('../results/blast/swissprot/%s.xml' % ngo)
        except:
            blasts = {}
        record = ngo + ':\n'
        try:
            with open('../data/ncbiSequence/%s.gb' % ngo, 'r') as f:
                ncbi = SeqIO.read(f, 'gb')
        except:
            continue
        try:
            with open('../data/uniprotSequence/%s.gb' % ngo, 'r') as f:
                uniprot = SeqIO.read(f, 'gb')
        except:
            uniprot = None
        if ncbi:
            record += '\tNCBI ID: %s\n' % ncbi.id
            record += '\tNCBI Name: %s\n' % ncbi.name
        if uniprot:
            record += '\tUniprot ID: %s\n' % uniprot.id
            record += '\tUniprot Name: %s\n' % uniprot.name
        if ngo in annotation_scores:
            record += '\tUniprot Annotation:\n'
            record += '\t\tRevision status: %s\n' % annotation_scores[ngo][0]
            record += '\t\tAnnotation score: %s\n' % annotation_scores[ngo][1]
            record += '\t\tAnnotation description: %s\n' % annotation_scores[ngo][2]
        if ngo in locations:
            record += '\tCellular Location:\n'
            for location in locations[ngo]:
                record += '\t\tLocation:\n'
                record += '\t\t\tlocation id: %s\n' % location[0]
                record += '\t\t\tlocation description: %s\n' % location[1]
        if ncbi:
            record += '\tDB Link:\n'
            for db in ncbi.dbxrefs[0].split(' '):
                db = db.split(':')
                record += '\t\t%s: %s\n' % (db[0], db[-1])
        if uniprot:
            for db in uniprot.dbxrefs[0].split(' '):
                db = db.split(':')
                record += '\t\t%s: %s\n' % (db[0], db[-1])
        if ncbi:
            with open('../data/ncbiGenome/sequence_full.gb', 'r') as f:
                genome = SeqIO.read(f, 'gb')
            for feature in genome.features:
                if feature.type == 'CDS':
                    if feature.qualifiers['locus_tag'][0] == ngo:
                        if feature.qualifiers['product'][0] in ref_ids:
                            record += '\tReference IDs: %s\n' % ', '.join(ref_ids[feature.qualifiers['product'][0]])
                        break
            record += '\tSequence: %s\n' % ncbi.seq
            record += '\tFeatures:\n'
            for feature in ncbi.features:
                record += '\t\tType: %s\n' % feature.type
                for qualifier in feature.qualifiers:
                    record += '\t\t\t%s: %s\n' % (qualifier, feature.qualifiers[qualifier][0])
                record += '\t\t\tlocation: %s\n' % str(feature.location).replace(':', '.').strip('order[]').strip('{}')
        if uniprot:
            for feature in uniprot.features:
                qualifiers = ''
                for qualifier in feature.qualifiers:
                    if qualifier == 'type':
                        tp = feature.qualifiers[qualifier][0][0].upper() + feature.qualifiers[qualifier][0][1:]
                    else:
                        qualifiers += '\t\t\t%s: %s\n' % (qualifier, feature.qualifiers[qualifier][0])
                qualifiers = '\t\t\tlocation: %s\n' % str(feature.location).replace(':', '.').strip('order').strip('[]{}')
                record += '\t\tType: %s\n' % tp + qualifiers
        if blasts:
            record += '\tBlast Results:\n'
            for alignment in blasts:
                record += '\tAlignment ID: %s\n' % alignment
                record += '\t\tdescription: %s\n' % blasts[alignment][0]
                record += '\t\tlength: %s\n' % blasts[alignment][1]
                record += '\t\te-value: %s\n' % blasts[alignment][2]
                record += '\t\tidentity: %s\n' % blasts[alignment][3]

        with open(output_file, 'a') as f:
            f.write(record + '\n')
