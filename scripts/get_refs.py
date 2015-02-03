#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez, SeqIO

Entrez.email = 'andremsantiago_@hotmail.com'
ltstart, ltend = 'NGO0971', 'NGO1212'
gi = '59800473'  # '59717368'

print('Fetching genome info...', end=' ', flush=True)
# request = Entrez.efetch(db="nucleotide", rettype='gb', retmode='text', id=gi)
request = open('./ncbiSeq/sequence_full.gb', 'r')
print('Done')
print('Reading genome...', end=' ', flush=True)
record = SeqIO.read(request, 'gb')
print('Done')
request.close()

with open('ids.txt', 'w') as f:
    f.write('')

db, names = 'pubmed', []

for feature in record.features:
    if feature.type != 'CDS':
        continue
    print(feature)
    locus_tag = feature.qualifiers['locus_tag'][0]
    product = feature.qualifiers['product'][0]
    if locus_tag < ltstart or locus_tag > ltend or product in names:
        continue
    names.append(product)
    term = 'Neisseria gonorrhoeae[Orgn] AND ' + product + '[Gene]'
    print('Searching number of references for gene', product, end='... ', flush=True)
    handle = Entrez.egquery(term=term)
    rec = Entrez.read(handle)
    handle.close()
    for row in rec['eGQueryResult']:
        if row['DbName'] == db:
            retmax = row['Count']
    print('Found', retmax)
    if retmax == '0':
        continue
    print('Fetching reference IDs', end='... ', flush=True)
    handle = Entrez.esearch(db=db, term=term, retmax=retmax)
    rec = Entrez.read(handle)
    idsList = ''
    for ID in rec['IdList']:
        idsList += ' ' + str(ID)
    with open('ids.txt', 'a') as f:
        f.write('Search term: ' + term + '\nIDs:' + idsList + '\n\n')
    print('Done')
