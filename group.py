#! /usr/bin/env python
#! /Users/user/anaconda3/bin/python3
#! /usr/bin/python3

#usage: ./group.py metadata.txt output.csv

import pandas as pd
import sys
from collections import Counter

df = pd.read_csv(sys.argv[1],sep="\t")
list_index=df.columns.tolist()
for i in range(len(list_index)):
    name = list_index[i]
    name = name + ".txt"
    name_list = df[list_index[i]].tolist()
    KOlist=[]
    for j in name_list:
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

combine.to_csv(sys.argv[2])
