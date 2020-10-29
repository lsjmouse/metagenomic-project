#! /usr/bin/env python
#! /Users/user/anaconda3/bin/python3
#! /usr/bin/python3

#usage: ./group.py metadata.txt output.csv

import pandas as pd
from collections import Counter
import argparse
usage = 'This programme combines all annotation result txt files into a csv table, and this programme needs to be the same directory of all txt files .A metadata txt is required for grouping principle.'
parser = argparse.ArgumentParser(description = usage)
parser.add_argument('-m', dest='metadata',metavar='METADATA',type=argparse.FileType('r'),required=True, help ='Metadata file for grouping principle')
parser.add_argument('-o',dest='outfile',metavar ='OUTFILE ',type = argparse.FileType('w') ,required=True, help = 'Output csv file')
args = parser.parse_args()

df = pd.read_csv(args.metadata,sep="\t")
list_index=df.columns.tolist()
for i in range(len(list_index)):
    name = list_index[i]
    name = name + ".txt"
    name_list = df[list_index[i]].tolist()
    name_cleaned_list = [x for x in name_list if str(x) != 'nan']
    KOlist=[]
    for j in name_cleaned_list:
        with open(j,"r") as fin:
            for line in fin:
                line=line.rstrip()
                if '\t' in line:
                    KOlist.append(line.split('\t')[1])

    a=Counter(KOlist)
    with open(name,"w") as fout:
        print('KO_ID\t{}'.format(list_index[i]), file = fout)
        for key in a:
            print('{}\t{}'.format(key,a[key]), file = fout)

firstname=list_index[0] + ".txt"
combine = pd.read_csv(firstname,sep="\t")
combine = combine.set_index("KO_ID")

for i in range(1,len(list_index)):
    name = list_index[i] + ".txt"
    df = pd.read_csv(name, sep="\t")
    df = df.set_index("KO_ID")
    combine = pd.merge(combine,df, how="outer",on=['KO_ID'])
    combine = combine.fillna(0)

combine.to_csv(args.outfile)
