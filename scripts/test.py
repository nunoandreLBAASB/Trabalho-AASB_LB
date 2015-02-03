# -*- coding: utf-8 -*-
#!/usr/bin/env python3


from Bio import SeqIO


# General info about sequence in analisys
def seq_info(record):
    seq_id = record.id   # protein ID
    seq_name = record.name  # protein Name
    seq_discr = record.description  # small description
    seq_size = len(record.seq)  # protein size
    seq_gb = record.seq  # aa protein sequence
    return seq_id, seq_name, seq_discr, seq_size, seq_gb


# Anotations about organism name and taxonomy
def annotation(record):
    seq_anota = record.annotations
    seq_organ = seq_anota['organism']  # organnism name
    seq_taxo = seq_anota['taxonomy']  # organism taxonomy
    return seq_organ, seq_taxo


#
def features(record):
    feats = []
    for i in range(len(record.features)):
        feats.append(record.features[i].type)  # , record.features[i].location)
    return feats


def qualifiers(record):
    qualifs = []
    for i in range(len(record.features)):
        if record.features[i].type == 'CDS':
            qualifs += (record.features[i].qualifiers['note'])
    return qualifs


# NAO CONSIGO POR A ABRIR TODOS OS FICHEIROS AUTOMATICAMENTE
def main():
    files = range(971, 1213)
    for file in files:
        try:
            with open('./SeqRecord/NGO' + ('0' if file < 1000 else '') + str(file) + '.gb', 'r') as request:
                record = SeqIO.read(request, 'gb')
                print(seq_info(record))
                # print(annotation(record))
                # print(features(record))
                # print(qualifiers(record))
        except IOError:
            pass

main()
