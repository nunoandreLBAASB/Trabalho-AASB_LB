#!/usr/bin/env python3


if __name__ == '__main__':
    # input file obtained from http://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz
    input_file = '../data/virulenceFactors/Neisseria_VFs_comparison.tsv'
    output_file = '../results/vfs.txt'
    species = 'N. gonorrhoeae FA 1090' # line[2], on line 28
    ltstart, ltend = 'NGO0971', 'NGO1212'
    tp = ''
    subtp = ''
    gene_name = ''

    with open(output_file, 'w') as f:
        pass  # clears file

    with open(input_file, 'r') as vfs:
        for line in vfs:
            line = line.rstrip().split('\t')
            counter = 0
            for entry in line:
                if len(entry) > 0:
                    counter += 1
            if counter == 1:
                tp = line[0]
                subtp = ''
                gene_name = ''
            else:
                if line[0]:
                    subtp = line[0]
                gene_name = line[1]
                locus_tag = line[2]
                if locus_tag >= ltstart and locus_tag <= ltend:
                    with open(output_file, 'a') as f:
                        f.write('%s\t%s\t%s\t%s\n' % (locus_tag, gene_name, tp, subtp))
