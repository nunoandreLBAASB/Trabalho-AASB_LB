#!/usr/bin/env python3


from Bio import SeqIO


def get_geneid(filename):
    geneid = ''
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    try:
        for db in record.features[-1].qualifiers['db_xref']:
            geneid = db.split(':')
            if geneid[0] == 'GeneID':
                geneid = geneid[1]
    except KeyError:
        geneid = '-'
    return geneid


def get_gene_accession(filename):
    accession = '-'
    with open(filename, 'r') as f:
        for line in f:
            try:
                start = line.index('Annotation: ') + len('Annotation: ')
                end = start + line[start:].index('.')
                accession = line[start:end]
                break
            except:
                continue
    return accession


def get_id(filename):
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    return record.id


def get_locustag(filename):
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    try:
        locus_tag = record.features[-1].qualifiers['locus_tag'][0]
    except KeyError:
        locus_tag = '-'
    return locus_tag


def get_gene_name(filename):
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    try:
        name = record.features[-1].qualifiers['gene'][0]
    except KeyError:
        name = '-'
    return name


def get_strand(filename, locus_tag):
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    strand = '-'
    for feature in record.features:
        if feature.type == 'gene' or feature.type == 'CDS':
            try:
                if feature.qualifiers['locus_tag'][0] == locus_tag:
                    strand = feature.location.strand
                    break
            except:
                pass
    return strand


def get_annotation_score(filename, locus_tag):
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip().split('\t')
            if line[0] == locus_tag:
                return line[1], line[2], line[3]
    return '-', '-', '-'


def get_protein_name(filename):
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    name = record.description
    return name


def get_number_aa(filename):
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    length = len(record.seq)
    return length


def get_location(filename, locus_tag):
    location = []
    found = False
    with open(filename, 'r') as f:
        for line in f:
            if len(line) > 1:
                line = line.rstrip().split('\t')
                if locus_tag == line[0] or found:
                    found = True
                    if line[1:]:
                        location.append(line[1:])
            elif found:
                break
    return location


def get_gos(filename):
    gos = []
    start = 0
    exist_gos = True
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    dbs = record.dbxrefs[0]
    while exist_gos:
        try:
            start += dbs[start:].index('GO:GO:') + len('GO:GO:')
            end = start + dbs[start:].index(' ')
            gos.append(dbs[start:end])
        except ValueError:
            exist_gos = False
    return gos


def get_ec(filename):
    ec = ''
    start = 0
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    dbs = record.dbxrefs[0]
    try:
        start = dbs[start:].index('EC:') + len('EC:')
        end = start + dbs[start:].index(' ')
        ec = dbs[start:end]
    except:
        ec = '-'
    return ec


def get_function_description(filename):
    with open(filename, 'r') as f:
        record = SeqIO.read(f, 'gb')
    try:
        note = record.features[-1].qualifiers['note'][0]
    except KeyError:
        note = '-'
    return note


def get_predicted_function(filename, locus_tag):
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.rstrip().split('\t')
                if line[0] == locus_tag:
                    return line[-1][line[-1].index(':') + 2:]
    except:
        pass
    return '-'


def get_vfs(filename, locus_tag):
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.rstrip().split('\t')
                if line[0] == locus_tag:
                    return line[1], line[2], line[3]
    except:
        pass
    return '-', '-', '-'


def get_directory_filenames(directory):
    import os

    for _, _, files in os.walk(directory):
        for filename in files:
            yield filename


if __name__ == '__main__':
    table = '../data/ProteinTable864_169534.tsv'
    genome = '../data/ncbiGenome/sequence_full.gb'
    gene = '../data/ncbiGene'
    ncbi = '../data/ncbiSequence'
    uniprot = '../data/uniprotSequence'
    score = '../results/score.txt'
    functions = '../results/functions.txt'
    vfs = '../results/vfs.txt'
    location_file = '../results/location.txt'
    output_file = '../results/table.tsv'
    columns = [ 'GeneID',
                'Accession',
                'Locus tag',
                'Gene name',
                'strand',
                'Uniprot ID',
                'Revision status',
                'Annotation score',
                'Annotation description',
                'ProteinID',
                'Protein Name',
                '# Aminoacid',
                'Virulence factor',
                'NCBI location ID',
                'GO\'s',
                'EC number',
                'Description',
                'Comments'
    ]

    with open(output_file, 'w') as f:
        f.write('\t'.join(columns) + '\n')

    for filename in get_directory_filenames(ncbi):
        print(filename.split('.')[0])
        geneid = get_geneid('%s/%s' % (ncbi, filename))
        gene_accession = '-'  # get_gene_accession('%s/%s.txt' % (gene, filename.split('.')[0]))
        locus_tag = get_locustag('%s/%s' % (ncbi, filename))
        gene_name = get_gene_name('%s/%s' % (ncbi, filename))
        strand = get_strand('%s' % genome, locus_tag)
        try:
            uniprot_id = get_id('%s/%s' % (uniprot, filename))
        except:
            uniprot_id = '-'
        annotation_score = get_annotation_score(score, locus_tag)
        protein_id = get_id('%s/%s' % (ncbi, filename))
        try:
            protein_name = get_protein_name('%s/%s' % (uniprot, filename))
            if protein_name == '.':
                protein_name = '-'
        except:
            protein_name = '-'
        try:
            aa = get_number_aa('%s/%s' % (uniprot, filename))
        except:
            aa = '-'
        vf = get_vfs(vfs, locus_tag)
        if gene_name == '-' and vf[0] != '-':
            gene_name = vf[0]
        vf = ';'.join(vf[1:])
        if vf == '-;-':
            vf = '-'
        location = get_location(location_file, locus_tag)
        for i in range(len(location)):
            location[i] = '%s (%s)' % (location[i][0], location[i][1])
        location = ';'.join(location)
        if not location:
            location = '-'
        try:
            gos = get_gos('%s/%s' % (uniprot, filename))
            gos = ';'.join(gos)
            if not gos:
                gos = '-'
        except:
            gos = '-'
        try:
            ec = get_ec('%s/%s' % (uniprot, filename))
        except:
            ec = '-'
        # if annotation_score[0] == '1':
        description = get_function_description('%s/%s' % (ncbi, filename))
        # # else:
        #     description = get_predicted_function(functions, locus_tag)
        #     if description != '-':
        #         description = 'Predicted from ' + description
        if annotation_score[0] == '1':
            comments = 'Reviewed by Swiss-Prot'
        else:
            comments = get_predicted_function(functions, locus_tag)
            if comments != '-':
                comments = 'Predicted from ' + comments


        line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (geneid, gene_accession, locus_tag, gene_name, strand, uniprot_id, annotation_score[0], annotation_score[1], annotation_score[2], protein_id, protein_name, str(aa), vf, location, gos, ec, description, comments)
        with open(output_file, 'a') as f:
            f.write(line + '\n')
