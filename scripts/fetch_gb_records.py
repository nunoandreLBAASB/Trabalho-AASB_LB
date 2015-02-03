#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez, SeqIO
Entrez.email = "andremsantiago_@hotmail.com"
gi, gis = '59800473', []
# lcs = {}


def import_seqrecord(gi):
    # fetches the full genome and converts to a SeqRecord
    request = Entrez.efetch(db="nucleotide", rettype='gbwithparts', retmode='text', id=gi)   # NCBI extract file (porque gbwithparts)
    record = SeqIO.read(request, 'gb')   # Read file from gb format to SeqRecord
    request.close()   # Close NCBI request
    return record


def read_seqrecord(filename, format):
    # Read a file
    with open(filename, 'r') as request:   # Open given file in read mode
        record = SeqIO.read(request, format)   # Use SeqIO function to read the file
    return record


genome = read_seqrecord('./ncbiSeq/sequence_full.gb', 'gb')   # Read file with all genome


# selects our relevant locus tags, between NGO0971 and NGO1212, and adds them to 'gis'
for feature in genome.features:
    if feature.type == 'CDS':
        try:
            locus_tag = feature.qualifiers['locus_tag'][0]   # Get the locus tags
            if locus_tag >= 'NGO0971' and locus_tag <= 'NGO1212':   # If between our target genes
                gis.append(feature.qualifiers['db_xref'][0][3:])   # Save the GI number of each gene
                # lcs[gis[-1]] = locus_tag
        except:
            pass

# fetches all the records according to the selected 'gis'
request = Entrez.efetch(db="nucleotide", rettype='gbwithparts', retmode='text', id=gis)   # NCBI request with select gis
records = SeqIO.parse(request, 'gb')   # Read request with SeqIO function from GB to individual SeqRecords in a list


# saves each individual SeqRecord to a 'gb' file, with its locus_tag as filename
for record in records:   # For each gi file
    locus_tag = record.features[-1].qualifiers['locus_tag'][0]   # From features.qualifiers select from CDS the locus_tag
    # with open('./ncbiSeqRecord/' + lcs[rec.annotations['gi']] + '.gb', 'w') as f:
    with open('./ncbiSeqRecord/%s.gb' % locus_tag, 'w') as f:
        SeqIO.write(record, f, 'gb')
