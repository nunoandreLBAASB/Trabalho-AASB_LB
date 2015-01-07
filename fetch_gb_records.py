#!/usr/bin/env python3

from Bio import Entrez, SeqIO
Entrez.email ="andremsantiago_@hotmail.com"
gi, gis = '59717368', []
lcs = {}

# fetches the full genome and converts to a SeqRecord
request = Entrez.efetch(db="nucleotide", rettype='gb', retmode='text', id=gi)
record = SeqIO.read(request, 'gb')
request.close()

# selects our relevant locus tags, between NGO0971 and NGO1212, and adds them to 'gis'
for feature in record.features:
	if feature.type == 'CDS':
		locus_tag = feature.qualifiers['locus_tag'][0]
		if locus_tag >= 'NGO0971' and locus_tag <= 'NGO1212':
			gis.append(feature.qualifiers['db_xref'][0][3:])
			lcs[gis[-1]] = locus_tag

# fetches all the records according to the selected 'gis'
request = Entrez.efetch(db="nucleotide", rettype='gb', retmode='text', id=gis)
records = SeqIO.parse(request, 'gb')

# saves each individual SeqRecord to a 'gb' file, with its locus_tag as filename
for rec in records:
	with open(lcs[rec.annotations['gi']] + '.gb', 'w') as f:
		SeqIO.write(rec, f, 'gb')
