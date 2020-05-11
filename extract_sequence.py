#! /usr/bin/env python
#! /Users/user/anaconda3/bin/python3
#! /usr/bin/python3

#this script is used to extract sequence from tblastx hit results.
#usage: ./extract_sequence.py input.fasta tblastx_outputfile output.fasta


import sys

with open(sys.argv[1],'r') as sequence_file, open(sys.argv[2],'r') as tblastx_file, open(sys.argv[3],'w') as fout:
    ID = []
    sequence = []
    for line in sequence_file:
        line = line.rstrip()
        if line.startswith('>'):
            id = line[1:]
            ID.append(id)
        else:
            sequence.append(line)
    sequence_dict = dict(zip(ID,sequence))

    tblastx_idset = set()
    for line in tblastx_file:
        id = line.split('\t')[0]
        tblastx_idset.add(id)

    for id in tblastx_idset:
        print('>{}\n{}'.format(id,sequence_dict[id]), file=fout)
