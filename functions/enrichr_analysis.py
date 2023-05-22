#!/usr/bin/python
from __future__ import division
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import seaborn as sns


def change_term(term):
    return term.split('_')[0]


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-inlist', dest='inlist', help='File containing gene-list')
    parser.add_argument('-out', dest='out', help='Folder for output')
    parser.add_argument('-q', dest='q', help='Run analysis with with or without enrichr query. True by default')
    args = parser.parse_args()

    #print(args.q)

    ##databasr api
    databases = ['TRANSFAC_and_JASPAR_PWMs','OMIM_Disease', 'OMIM_Expanded', 'KEGG_2016', 'BioCarta_2016','WikiPathways_2016','Panther_2016','GO_Biological_process_2018','GO_Cellular_Component_2018','GO_Molecular_Function_2018','Human_Phenotype_Ontology','MGI_Mammalian_Phenotype_2017', 'Jensen_DISEASES']
    #args.out = '../data/%s'+args.out
    lModules = []
    lFailed = []
    lOutfiles = []

    if args.out not in os.listdir('.'):
        if args.out not in os.listdir('.'):
            cmd = 'mkdir %s' %args.out
            os.system(cmd)

    #Extract Enrichr analysis
    for db in databases:
        out = db+'_'
        fOut = '%s/%s.txt' %(args.out,out)
        lOutfiles.append(fOut)
        if args.q == '1':
            cmd = 'python3 ./enrichr-api/query_enrichr_py3.py %s %s %s %s/%s' %(args.inlist, out, db, args.out, out)
            print(cmd)
            os.system(cmd)

    lOut_dbs = []
    lGene_members = []
    
    #Plot significant terms
    for out_db in lOutfiles:
        iFileLength = file_len(out_db)
        df = pd.read_csv(out_db, sep='\t', index_col=None)
        df = df.sort_values(by='Adjusted P-value')
        len_df_10 = len(df.loc[df['Adjusted P-value']<= 0.05])

        if len_df_10 >= 20:
            df_20 = df.head(20)
        else:
            df_20 = df.head(len_df_10)
        
        if len_df_10>0:
            df_10 = df.head(len_df_10)
            df_10['-log(padj)'] = df_10['Adjusted P-value'].map(lambda a: -np.log(a))
            df_10['Genes involved (%)'] = df_10['Overlap'].map(lambda a: 100*(int(a.split('/')[0])/int(a.split('/')[1])))
            
            if 'KEGG' in out_db:
                df_10['Term'] = df_10['Term'].apply(change_term)
            #plt.subplots(figsize=(14,8))
            plt.subplots(figsize=(14,len_df_10))
            ax = sns.barplot(x=df_10['Genes involved (%)'], y=df_10['Term'], orient='h', color='#998ec3')
            ax.set_ylabel('')
            #ax.set_axis_bgcolor('white')
            ax.set_fc('white')
            ax.set_xlabel('Genes involved (%)', size=18)
            plt.tick_params(axis='both', labelsize=16)
            ax2 = ax.twiny()
            ax2.set_xlabel('-log(padj)', size=18)
            ax2.yaxis.grid(b=None)
            ax2.grid(False)
            lineColor = '#f1a340'
            f = lambda x: l2[0] + (x - l[0]) / (l[1] - l[0]) * (l2[1] - l2[0])
            ax2.plot(sorted(df_10['-log(padj)']), range(len_df_10-1, -1, -1), marker='o', color=lineColor, linewidth=4,
                               markersize=10, markerfacecolor=lineColor, markeredgecolor=lineColor)

            plt.tick_params(axis='both', labelsize=16)
            fSaveplot = out_db.split('.txt')[0]
            plt.savefig('%s.pdf' %(fSaveplot), bbox_inches='tight')
        else:
            print('No significant terms', out_db)
    if lFailed:
        print('Failed: ', lFailed)

