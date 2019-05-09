import sys, csv, pandas as pd, json, glob
from pandas.io.json import json_normalize
#https://www.quora.com/How-does-one-convert-fasta-to-CSV-using-python

appended_data = []

for n, jsons in enumerate(glob.glob('{}/*.json'.format(sys.argv[1]))):

    print('\t Found json file {}: {}...'.format(n+1, jsons))

    with open(jsons) as json_file:  
        data = json.load(json_file)
        normalised = json_normalize(data)
        transposed = normalised.T
        transposed = transposed.reset_index()
        transposed['ORF_no'] = jsons
        transposed.ORF_no = transposed.ORF_no.replace(r".*\/","")
        transposed.columns = ['gene','results', 'ORF_no']
        transposed.gene = transposed.gene.replace(".*results\.","", regex=True)
        transposed.ORF_no = transposed.ORF_no.replace(".*_(.*)\..*","\\1", regex=True)
        appended_data.append(transposed)
        print(transposed)

appended_data = pd.concat(appended_data, axis=0)

#appended_data = appended_data.loc[appended_data.iloc[:,1].str.contains(r'^Escherichia', na=False)]

ORF_filter = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 44]

appended_data.ORF_no = appended_data.ORF_no.astype(int)

appended_data = appended_data[appended_data['ORF_no'].isin(ORF_filter)]

appended_data.to_csv("{}_filtered.csv".format(sys.argv[2]))