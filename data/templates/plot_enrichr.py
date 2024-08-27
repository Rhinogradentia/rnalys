#!/usr/bin/env python
# coding: utf-8


#script for enrichr plot

#input enrichr table
#output script generating plot

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import re
import os

def plot_enrichr(df, name_out=None, outfolder=None, title=None, pvalue_type='padj', figsize=(8,3), color=None, nterms=10, fontsize=10, exclude_go=False, sort_by='pvalue'):
    
    
    if color == None:
        color = 'yellow'
    
    nterms = min(nterms, df.shape[0])
    df = df.head(nterms)

    color = re.sub(r'\d+', '', color)

    plt.subplots(figsize=figsize)
    ax = sns.barplot(x=df['Genes involved (%)'], y=df['Term'], orient='h', color=color)
    ax.set_ylabel('')
    #ax.set_axis_bgcolor('white')
    ax.set_xlabel('Genes involved (%)', size=fontsize)
    plt.tick_params(axis='both', labelsize=fontsize)
    ax2 = ax.twiny()
    ax2.set_xlabel('-log(%s)' %pvalue_type, size=fontsize)
    ax2.yaxis.grid(b=None)
    ax2.grid(False)
    lineColor = '#fc8d59'
    f = lambda x: l2[0] + (x - l[0]) / (l[1] - l[0]) * (l2[1] - l2[0])
    #ax2.plot(df['-log(%s)'%pvalue_type], range(nterms-1, -1, -1), marker='o', color=lineColor, linewidth=2,
    #                   markersize=5, markerfacecolor=lineColor, markeredgecolor=lineColor)

    ax2.plot(df['-log(%s)'%pvalue_type].iloc[::-1], range(nterms-1, -1, -1), marker='o', color=lineColor, linewidth=3, markersize=7, markerfacecolor=lineColor, markeredgecolor=lineColor)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    plt.axvline(x=-np.log(0.05), color='red')

    #print(df)
    if title:
        plt.title(title)

    plt.tick_params(axis='both', labelsize=fontsize)


def sort_filter_df(df, sort_by='padj', only_significant=True, pvalue_type='padj', filter_go=False):
    '''
    
    For sorting and filtering enricher dataframes
    
    '''
    
    #df = pd.read_csv(out_db, sep='\t', index_col=None)
    if pvalue_type == 'padj':
        pvalue_type_indf = 'Adjusted P-value'
        
    elif pvalue_type == 'pval':
        pvalue_type_indf = 'P-value'
        
    #print(df)    
    if only_significant:
        df =  df.loc[df[pvalue_type_indf]<=0.05,]
    
    if sort_by is not None:
        df = df.sort_values(by=pvalue_type_indf)
    
    df.loc[:,'-log(%s)' %pvalue_type] = df.loc[:, pvalue_type_indf].map(lambda a: -np.log(a))
    df.loc[:,'Genes involved (%)'] = df.loc[:,'Overlap'].map(lambda a: 100*(int(a.split('/')[0])/int(a.split('/')[1])))
    
    if sort_by == 'percgenesinvolved':
        df = df.sort_values(by='Genes involved (%)', ascending=False)

    elif sort_by== 'genesinvolved':
        df['genesinvolved'] = df.loc[:, 'Overlap'].map(lambda a: int(a.split('/')[0]))
        df = df.sort_values(by='genesinvolved', ascending=False)
        
    elif sort_by == 'padj':
        df = df.sort_values(by='Adjusted P-value', ascending=True)
        
    elif sort_by == 'Combined Score':
        df = df.sort_values(by='Combined Score', ascending=False)
        
    else:
        #df = df.sort_values(by=sort_by, ascending=False)
        #try:
        df = df.sort_values(by=sort_by, ascending=False)
        #except
        #    print('Enter a correct sort term')
        #    break
        
    if filter_go == True:
        df.Term = [x.split('(GO:')[0] for x in df.Term if '(GO' in x]
    
    return df


#os.chdir('')

folder = '../data/generated/enrichr/'
for file_ in os.listdir(folder):
    if '_fileid_' in file_:
        infile = os.path.join(folder, file_)
        df = pd.read_csv(infile, sep='\t')
        title = infile.split('/')[-1]
        plot_enrichr(sort_filter_df(df), title=title, color='yellow', name_out=title+'.pdf')






