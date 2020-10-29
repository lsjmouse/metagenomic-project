#! /usr/bin/env python
#! /Users/user/anaconda3/bin/python3
#! /usr/bin/python3

import argparse
usage = 'This programme adds 50 N to the unassembled sequence files from Pear'
parser = argparse.ArgumentParser(description = usage)
parser.add_argument('-f', dest='forward',metavar='Forward Sequence',type=argparse.FileType('r'),required=True, help ='Input forward fastq')
parser.add_argument('-r', dest='reverse',metavar='Reverse Sequence',type=argparse.FileType('r'),required=True, help ='Input reverse fastq')
parser.add_argument('-o',dest='outfile',metavar ='OUTFILE ',type = argparse.FileType('w') ,required=True, help = 'Output file')
args = parser.parse_args()


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'}
    return ''.join([complement[base] for base in dna[::-1]])


counter1 = 0
forwardid = []
forwardseq = []
forwardwholeid = []
for line in args.forward:
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
for line in args.reverse:
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
    print('>{}\n{}{}{}'.format(iddict[id][1:],forwarddict[id],insert,reversedict[id]), file=args.outfile ,sep='')
