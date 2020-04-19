#! /Users/user/anaconda3/bin/python3
#! /usr/bin/python3

#descrption: use for combine unassembled sequence from Pear with adding 50 'N'
#usage: ./pear_add_n.py forward.fastq reverse.fastq output.fasta


import sys

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'}
    return ''.join([complement[base] for base in dna[::-1]])

with open(sys.argv[1],'r') as forward, open(sys.argv[2],'r') as reverse, open(sys.argv[3],'w') as fout:
    counter1 = 0
    forwardid = []
    forwardseq = []
    forwardwholeid = []
    for line in forward:
        line = line .rstrip()
        if counter1%4 == 0:
            forwardwholeid.append(line)
            id = line.split()[0]
            forwardid.append(id)
            counter1 += 1
        elif counter1%4 == 1:
            forwardseq.append(line)
            counter1 += 1
        else:
            counter1 += 1
    forwarddict=dict(zip(forwardid,forwardseq))
    iddict=dict(zip(forwardid,forwardwholeid))

    counter2 = 0
    reverseid = []
    reverseseq = []
    for line in reverse:
        line = line .rstrip()
        if counter2%4 == 0:
            id=line.split()[0]
            reverseid.append(id)
            counter2 += 1
        elif counter2%4 == 1:
            reverse_complement_seq = reverse_complement(line)
            reverseseq.append(reverse_complement_seq)
            counter2 += 1
        else:
            counter2 += 1
    reversedict=dict(zip(reverseid,reverseseq))


    insert = 50*'N'
    for id in forwardid:
        print('>{}\n{}{}{}'.format(iddict[id][1:],forwarddict[id],insert,reversedict[id]), file=fout,sep='')
