#! /usr/bin/env python
#! /Users/user/anaconda3/bin/python3
#! /usr/bin/python3
import sys
def reverse(sequence):
    sequence=sequence.upper()
    old_chars='ACGT'
    new_chars='TGCA'
    tab=str.maketrans(old_chars,new_chars)
    return sequence.translate(tab)[::-1]
def transcribe(sequence):
    sequence=sequence.upper()
    for ch in sequence:
        rna_seq = sequence.replace('T', 'U')
        return(rna_seq)
def translate(sequence, position):
    codon = {"AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N",
            "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T",
            "AGA":"R", "AGC":"S", "AGG":"R", "AGU":"S",
            "AUA":"I", "AUC":"I", "AUG":"M", "AUU":"I",
            "CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H",
            "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P",
            "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R",
            "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L",
            "GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D",
            "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A",
            "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G",
            "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V",
            "UAA":"*", "UAC":"Y", "UAG":"*", "UAU":"T",
            "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S",
            "UGA":"*", "UGC":"C", "UGG":"W", "UGU":"C",
            "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F"}
    protein_seq=''
    for n in range(int(position), len(sequence), 3):
        if sequence[n:n+3] in codon:
            protein_seq += codon[sequence[n:n+3]]
    return protein_seq
ids, sequences = [], []
memoery=''
with open(sys.argv[1]) as fin,open(sys.argv[2],'w') as fout:
    for line in fin:
        line = line.rstrip()
        if line[0] == '>':
            if memoery=='':
                id = line.split(' ')[0][1:]
                #id = line.split('\t')[0][1:]
                ids.append(id)
            else:
                sequences.append(memoery)
                memoery=''
                id = line.split('\t')[0][1:]
                ids.append(id)
        else:
            memoery += line
    sequences.append(line)
    dict_protein=dict(zip(ids,sequences))
    for key, value in dict_protein.items():
        rna = transcribe(value)
        protein_1 = translate(rna,0)
        protein_2 = translate(rna,1)
        protein_3 = translate(rna,2)
        reverse_dna=reverse(value)
        reverse_rna=transcribe(reverse_dna)
        reverse_protein_1=translate(reverse_rna,0)
        reverse_protein_2=translate(reverse_rna,1)
        reverse_protein_3=translate(reverse_rna,2)
        print('>{}:1\n{}\n>{}:2\n{}\n>{}:3\n{}\n>{}:4\n{}\n>{}:5\n{}\n>{}:6\n{}'.format(
        key,protein_1,key,protein_2,key,protein_3,key,reverse_protein_1,key,reverse_protein_2,
        key,reverse_protein_3), file=fout)
