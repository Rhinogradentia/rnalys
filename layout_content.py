
import os
import json
import math
import flask
import dash
import dash_auth
import dash_table as dt
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import dash_bio as dashbio
import pyfiglet
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash import dcc
from dash import html
from dash import Input, Output, State
from dash import dash_table
from collections import Counter
from dash.exceptions import PreventUpdate
from sklearn.decomposition import PCA
from demos import dash_reusable_components as drc



#Genesymbol
df_symbol_file = './data/ensembl_symbol.csv'
df_symbol = pd.read_csv(df_symbol_file, sep='\t')
df_symbol.index = df_symbol['ensembl_gene_id']
df_symbol = df_symbol.loc[~df_symbol.index.duplicated(keep='first')]
df_symbol_sym_ind = df_symbol.copy()
df_symbol_sym_ind.index = df_symbol_sym_ind['hgnc_symbol']
df_symbol['hgnc_symbol'].fillna(df_symbol['ensembl_gene_id'], inplace=True)
dTranslate = dict(df_symbol['hgnc_symbol'])

dups = [k for k, v in Counter(df_symbol['hgnc_symbol'].tolist()).items() if v > 1]
for k, v in dTranslate.items():
    if v in dups:
        dTranslate[k] = k

df_symbol_sym_ind.index = df_symbol_sym_ind.index.fillna('NaN')
hgnc_dropdown = list(df_symbol_sym_ind.index)
df_counts_combined_file = '/home/cfrisk/Dropbox/dash/data/df_counts_combined.tab'
df_counts_combined = pd.read_csv(df_counts_combined_file, sep='\t', index_col=0)

df_meta_combined_file = '/home/cfrisk/Dropbox/dash/data/df_meta_v3.tab'
df_meta_combined = pd.read_csv(df_meta_combined_file, sep='\t', index_col=0)
available_tissues = df_meta_combined['tissue'].unique()

start_genes = ['LEP', 'ADIPOQ', 'RARRES2', 'IL6', 'CCL2', 'SERPINE1', 'RBP4', 'NAMPT', 'ITLN1', 'GRN']

#for volcano init
columns = ['Ensembl', 'hgnc', 'baseMean', 'log2FoldChange', 'pvalue', 'padj']
data = [['ens', 'hgnc1', 10, 4, 0.005, 0.9995]]
df = pd.DataFrame(data, columns=columns)


#Load datasets
fDatasets = 'data/datasets.txt'

if os.path.isfile('data/datasets.txt'):
    lPapers = []
    with open(fDatasets, 'r') as f:
        for line in f:
            sPaper = line.split()[0]
            lPapers.append(sPaper)
    lPapers.append('New')
else:
    os.system('touch data/datasets.txt')

layout_index = html.Div([
    html.Div([
       html.P('Setup')
    ]),
    html.Div([
        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Count File')
            ]),
            style={
                'width': '40%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            # Allow multiple files to be uploaded
            multiple=False
        ),
        dcc.Store(id='store_counts'),
        html.Div(id='output-data-upload'),

        dcc.Upload(
            id='upload-data_info',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Info File')
            ]),
            style={
                'width': '40%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            # Allow multiple files to be uploaded
            multiple=False
        ),
        dcc.Store(id='store_info'),
        html.Div(id='output-data-upload-info'),

    ]),

    dcc.Link('DE + PCA + Enrichr', href='/page-1'),
    html.Br(),
    #dcc.Link('WGCNA + Enrichr', href='/page-2'),
])

layout_page1 = html.Div(
    children=[
        # Error Message
        html.Div(id="error-message"),
        # Top Banner
        html.Div(
            className="study-browser-banner row",
            children=[
                html.H2(className="h2-title", children="RNA-analysis"),
            ],
        ),
        html.Div([
            html.Div([
                html.P('Load dataset:',
                       style={'width': '30%', 'display': 'flex', 'verticalAlign': "middle", 'padding': '2px'}),
                html.Div([
                    dcc.Dropdown(id='paper',
                                 options=[{'label': j, 'value': j} for j in lPapers],
                                 value='New')
                ], style={'width': '30%', 'display': 'inline-block'}),
            ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '50%'}),

            html.Div([
                dbc.Input(
                    id='save_dataset_name',
                    type='text',
                    placeholder='name dataset'
                ),
                dbc.Button('Save', id='btn-save_dataset', n_clicks=0),
                html.Div(id='name_dataset_save'),
            ], style={'width': '30%', 'display': 'inline-block', 'padding':'5px'}),


        ]),

        html.P(id='dataset_name_placeholder'),
        html.Div([
            html.Div([
                html.Div([
                    html.H3('Select samples', style={'textAlign': 'left', 'padding': '5px'}),

                    html.Div([
                        html.P('Type:', style={'width': '10%', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),

                        html.Div([
                            dcc.Dropdown(
                                id='type_dropdown',
                                options=[{'label': j, 'value': j} for j in df_meta_combined['type'].unique().tolist() + ['HF']],
                                value=[''],
                                multi=True,
                            )
                        ], style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),

                    html.Div([
                        html.P('Tissue:',style={'width': '10%', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),
                        html.Div([
                            dcc.Dropdown(
                            id='tissue_dropdown',
                            options=[
                                {'label': 'Left Ventricle', 'value': 'LV'},
                                {'label': 'Right Ventricle', 'value': 'RV'},
                                {'label': 'Right atrial', 'value': 'RAA'},
                                {'label': 'Muscle', 'value': 'musc'},
                                {'label': 'EpiFat', 'value': 'EF'},
                                {'label': 'SubFat', 'value': 'SF'},
                                {'label': 'LV+RV', 'value': 'LVRV'},
                                {'label': 'EF+SF', 'value': 'EFSF'}
                            ],
                            value=[''],
                            multi=True)
                        ], style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),
                    #
                    html.Div([
                        html.P('Batch:', style={'width': '10%', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),
                        html.Div([
                            dcc.Dropdown(
                            id='batch_dropdown',
                            options=[{'label': j, 'value': j} for j in df_meta_combined['SeqTag'].unique()],
                            multi=True)
                        ], style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),

                    html.Div([
                        html.Div(id='exclude_list', style={'display': 'none'}),
                        html.P('Exclude:', style={'margin-right': '10px'}),
                        html.Div([
                            dcc.Dropdown(
                                id='exclude',
                                options=[{'label': j, 'value': j} for j in df_counts_combined.columns],
                                multi=True,
                                placeholder='Select Samples to exclude',

                            )
                        ],style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),

                    #html.Div(id='read_table'),
                    html.Div(id='selected_data', style={'display': 'none'}),
                    html.Div(id='intermediate-table', style={'display': 'none'}),

                ], className='row', style={'width': '40%','display': 'inline-block', 'margin':'10px'}),

                html.Div([
                    html.H3('Transformation & Normalization', style={'textAlign': 'left', 'padding': '10px'}),

                    html.Div([
                           html.Div([
                               html.P('Select transf/norm:',
                                      style={'width': '30%', 'display': 'flex', 'verticalAlign': "middle",
                                             'padding': '2px'}),
                               html.Div([
                                   dcc.Dropdown(
                                       id='transformation',
                                       options=[{'label': j, 'value': j} for j in ['Sizefactor normalization', 'vsd','rlog','quantilenorm_log2']],
                                       multi=False,
                                       placeholder='Select Samples to exclude',)
                               ], style={'display': 'inline-block', 'width': '100%', 'verticalAlign': "middle"})
                           ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),
                            #], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),

                        html.Div([
                            html.P('Remove Confounding:', style={'width': '30%', 'display': 'flex', 'verticalAlign':"middle", 'padding':'2px'}),
                            html.Div([
                                dcc.Dropdown(
                                    id='rm_confounding',
                                    options=[{'label': j, 'value': j} for j in df_meta_combined.columns],
                                    multi=False,
                                    value=None)
                            ], style={'display': 'inline-block', 'width': '100%', 'verticalAlign': "middle"})
                        ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),



                        #)

                        html.Div([
                            dcc.Checklist(
                                id='force_run',
                                options=[{'label': ' Force re-run', 'value': 'force_run'}],
                                labelClassName='force_run'
                            )
                        ]),
                        html.Br(),

                        #dbc.Button("Submit", id='btn-selected_data_submit', color="primary", className="me-1",
                        #           style={'display': 'block', 'margin': 'auto', 'width':'30%'}),
                        html.Div(
                            [
                                html.Div([
                                    dbc.Button('Submit', id='btn-selected_data_submit', n_clicks=0),
                                ], style={'width':'40%'}),
                                dbc.Spinner(html.Div(id="submit-done"),size='md', delay_hide=1, delay_show =1, show_initially=False),
                            ], style={'width':'60%'}
                        )
                    ]),
                ], className='row', style={'width': '50%', 'display': 'inline-block'})

            ], className='row', style={'display': 'flex', 'margin': 'auto', 'border': '1px solid #C6CCD5', 'padding': '20px', 'border-radius': '5px','justify-content': 'space-between', 'margin-top': 5, 'padding': 5}),
        #]),


            html.Div([
                dmc.Alert(
                    "-",
                    id="alert_selection",
                    color='info',
                    #is_open=True,
                    withCloseButton=True,
                    # n_clicks=0,
                ),
            ], style={'textAlign': 'left', 'width': '30%', 'margin-top': 5}),

            html.Div([
                dmc.Alert(
                    "-",
                    id="alert_main_table",
                    color='info',
                    withCloseButton=True
                    #n_clicks=0,
                ),
            ], style={'textAlign': 'left', 'width':'30%', 'margin-top': 5}, id='alert_main_toggle'),

        ], className='row', style={'display': 'flex', 'justify-content': 'space-between', 'margin-top': 5, 'padding': 5}),
        html.Br(),
        html.Div([
            html.Div([
                html.P('Select gene'),

                html.Div(id='Select_gene_input'),
                html.Div([
                    html.Div([html.P('Ensembl:')], style={'width': '5%', 'display': 'inline-block'}),
                    html.Div([dcc.Dropdown(id='input-2',
                                 options=[{'label': j, 'value': j} for j in df_counts_combined.index], multi=True,
                                 value=['ENSG00000232810'])],
                                 style={'width': '70%', 'display': 'inline-block', 'verticalAlign': "middle"})
                ]),

                html.Div([
                    html.Div([html.P('hgnc:')], style={'width': '5%', 'display': 'inline-block'}),
                    html.Div([dcc.Dropdown(id='input-2_hgnc', options=[{'label': j, 'value': j} for j in hgnc_dropdown],
                                 multi=True, value=start_genes)],style={'width': '70%', 'display': 'inline-block', 'verticalAlign':"middle"})
                ]),
            ], style={'display': 'inline-block', 'width':'70%', 'verticalAlign':"middle"}),




            html.Div([

                html.Div([
                    html.P('Display options'),

                    html.Div([
                        dcc.RadioItems(id='radio_symbol',
                                       options=[{'label': j, 'value': j} for j in ['Ensembl', 'Symbol']],
                                       value='Ensembl', labelStyle={'display': 'inline-block'}),
                    ]),
                    html.Div([
                        dcc.RadioItems(id='radio-grouping',
                                       options=[
                                           {'label': 'Tissue', 'value': 'tissue'},
                                           {'label': 'Type', 'value': 'type'},
                                           {'label': 'Batch', 'value': 'SeqTag'},
                                       ],
                                       value='tissue', labelStyle={'display': 'inline-block'}),
                    ]),

                    dcc.Checklist(id='full_text',
                                  options=[
                                      {'label': 'Full text', 'value': 'fulltext'},
                                  ], labelClassName='Full text'),

                ], style={'verticalAlign':"middle"}

                ),

            ], style={'margin-top': 10, 'border':
                '1px solid #C6CCD5', 'padding': 5,
                      'border-radius': '5px', 'display':'inline-block', 'width': '30%', 'verticalAlign':"middle"}
            ),
        ], className='row', style={'width': '100%','display': 'inline-block', 'verticalAlign':"middle"}),

        html.Div(id='log2_table', style={'display': 'none'}),

        #TABS
        html.Div([
                dbc.Tabs(
                    [
                        dbc.Tab(label="Gene Expression", tab_id="tab-1"),
                        dbc.Tab(label="Human protein atlas", tab_id="tab-2"),
                    ],
                    id="tabs",
                    active_tab="tab-1",
                ),
                html.Div(id="content"),
        ]),

        #dcc.Graph(id='indicator-graphic2'),
        #dcc.Graph(id='hpa_graph'),
        #dash_table.DataTable(id='gene_table',
        #                     columns=[
        #                         {"name": i, "id": i} for i in sorted(df_symbol.columns)
        #                     ]),

        html.Br(),
        html.Div([html.H2('PCA analysis')], style={"textAlign": "left"}),


        dbc.Input(id='input-gene', type='text', placeholder='Insert Gene (Ensembl)'),
        dbc.Input(id='num_genes', type='text', placeholder='# genes'),

        html.Div([
            dcc.Dropdown(
                id='meta_dropdown',
                options=[{'label': j, 'value': j} for j in df_meta_combined.columns],
                value='tissue'
            ),
        ], style={'width':'60%'}),

        dcc.RadioItems(id='radio_coloring', options=[{'label': j, 'value': j} for j in ['tissue', 'type', 'batch',
                                                                                        'Sex', 'gene-level']],
                       value='tissue', labelStyle={'display': 'inline-block'}),

        dcc.RadioItems(id='radio-text', options=[{'label': j, 'value': j} for j in ['None', 'Text']], value='None',
                       labelStyle={'display': 'inline-block'}, labelClassName='Scatter Text'),

        dcc.Checklist(
            id='biplot_radio',
            options=[{'label': 'Biplot', 'value': 'biplot'}],
        ),

        dcc.RadioItems(id='biplot_text_radio'),


        html.Div(children=[

                          dcc.Graph(id='pca_and_barplot', style={'display': 'inline-block', 'width':'50%'}),
                          dcc.Graph(id='barplot', style={'display': 'inline-block', 'width':'50%'})
            ]),

        html.Div(children=[

            dcc.Graph(id='pca-graph', style={'display': 'inline-block', 'width': '50%'}),
            dcc.Graph(id='pca_comp_explained', style={'display': 'inline-block', 'width': '50%'})
        ]),

        html.Div(id='clicked'),

        html.Div([
            html.Div([dcc.Graph(id='clustergram')]),
            html.Div([dcc.Graph(id='pca_correlation')]),
        ]),





        html.Div([
            html.H3('DE analysis'),

            html.Div([
                html.P('Settings'),



                html.Div([
                    #dcc.RadioItems(id='radio-de', options=[{'label': j, 'value': j} for j in ['DESeq2', 'edgeR', 'limma']],
                    #       value='DESeq2',
                    #       labelStyle={'display': 'inline-block'}),

                    html.Div([
                        html.P('Program:',
                               style={'width': '10%', 'display': 'flex', 'verticalAlign': "middle", 'padding': '5px'}),

                        html.Div([
                            dcc.Dropdown(
                                id='program',
                                options=[{'label': j, 'value': j} for j in ['DESeq2', 'edgeR', 'limma']],
                                value = [''],
                                multi = False,
                            )
                        ], style = {'display': 'inline-block', 'width': '20%', 'verticalAlign': "middle"})
                    ], style = {'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),




                    html.Div([
                        html.P('Rowsum:',
                               style={'width': '10%', 'display': 'flex', 'verticalAlign': "middle", 'padding': '5px'}),
                        html.Div([
                            dbc.Input(id='rowsum', type='text', value=10,
                                      placeholder='Select rowsum')
                        ], style={'display': 'inline-block', 'width': '20%', 'verticalAlign': "middle"})
                    ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),

                    html.Div([
                        html.P('Desing formula:',
                               style={'width': '10%', 'display': 'flex', 'verticalAlign': "middle", 'padding': '5px'}),
                        html.Div([
                            dbc.Input(id='design', type='text', value='~',
                                      placeholder='')
                        ], style={'display': 'inline-block', 'width': '20%', 'verticalAlign': "middle"})
                    ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),

                    html.Div([
                        html.P('Reference:',
                               style={'width': '10%', 'display': 'flex', 'verticalAlign': "middle", 'padding': '5px'}),
                        html.Div([
                            dbc.Input(id='reference', type='text',
                                      placeholder='Enter reference for DE')
                        ], style={'display': 'inline-block', 'width': '20%', 'verticalAlign': "middle"})
                    ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),

                ]),

                html.Div([
                    html.Div([
                        html.P('Display, X, Y',
                               style={'width': '15%', 'display': 'flex', 'verticalAlign': "middle", 'padding': '5px'}),
                        html.Div([
                            dcc.Dropdown(
                                id='volcano_xaxis',
                                options=[{'label': j, 'value': j} for j in ['log2FoldChange', 'baseMean']],
                                value=['log2FoldChange'],
                                multi=False,
                                placeholder='x-axis'
                            )], style={'width': '20%', 'display': 'inline-block', 'verticalAlign': "middle"}),

                        html.Div([dcc.Dropdown(
                            id='volcano_yaxis',
                            options=[{'label': j, 'value': j} for j in ['log2FoldChange', '-log10(p)', 'baseMean']],
                            value=['-log10(p)'],
                            multi=False,
                            placeholder='y-axis'
                        )], style={'width': '20%', 'display': 'inline-block', 'verticalAlign': "middle"}),
                    ], style={'width': '50%', 'display': 'flex', 'verticalAlign': "middle"}),
                ]),
                html.Div(
                [
                    #html.Div([
                    dbc.Button('Run analysis', id='btn-DE', n_clicks=None),
                    #], style={'width': '40%', 'padding': '5px'}),

                    dbc.Button('Export table', id='btn_export'),
                    dbc.Button('Export plot', id='btn_export_plot'),
                    dbc.Spinner(html.Div(id="submit-de-done"), size='md', delay_hide=1, delay_show=1,
                                show_initially=False),
                ], style={'width': '100%', 'padding':'5px'}
                ),
                #html.Div([

                #], style={'width': '100%', 'padding':'5px'})
            ], style={'margin-top': 10, 'border':
                                    '1px solid #C6CCD5', 'padding': 5,
                                    'border-radius': '5px'}),


            html.Br(),



            html.P(id='export_placeholder'),
            html.Div(id='intermediate-DEtable', style={'display': 'none'}),
            html.Div(id='temp', style={'display': 'none'}),
            html.Div(id='pvalue', style={'display': 'none'}),
        ]),

        html.Div(id='export_plot_clicked'),

        html.Div([
            'Effect sizes',
            dcc.RangeSlider(
                id='volcanoplot-input',
                min=-8,
                max=8,
                step=0.05,
                marks={
                    i: {'label': str(i)} for i in range(-8, 8)
                },
                value=[-1, 1]
            ),
            html.Br(),
            html.Div(
                dcc.Graph(
                    id='volcanoplot',
                    figure=dashbio.VolcanoPlot(effect_size='log2FoldChange', p='padj', gene='Ensembl',
                                               logp=True, snp='Ensembl', xlabel='log2FoldChange',
                                               genomewideline_value=2.5, dataframe=df
                    )
                )
            )
        ]),
        html.Div(id='number_of_degenes'),
        html.Div([


            html.Div([
                html.P('Filter on significance:',
                       style={'width': '15%', 'display': 'flex', 'verticalAlign': "middle", 'padding': '5px'}),
                html.Div([
                    dbc.Input(id='toggle_sig', type='text', value='0.05',
                              placeholder='Select p-value')
                ], style={'display': 'inline-block', 'width': '20%', 'verticalAlign': "middle"})
            ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),

            html.Div([
                html.P('Basemean',
                           style={'width': '15%', 'display': 'flex', 'verticalAlign': "middle", 'padding': '5px'}),
                html.Div([
                    dbc.Input(id='toggle_basemean', type='text', value='0',
                              placeholder='basemean')
                ], style={'display': 'inline-block', 'width': '20%', 'verticalAlign': "middle"})
            ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),

            html.Div([
                dbc.Button(id='sig_submit', n_clicks=0, children='Submit'),
            ], style={'display': 'inline-block', 'verticalAlign':"middle"}),


            html.Div([dbc.Button('Add DE set', id='add_de_table')],
                     style={'display': 'inline-block', 'verticalAlign':"middle"}),
            html.Div([
                dbc.Button('clear DE set', id='clear_de_table')],
                    style={'display': 'inline-block', 'verticalAlign':"middle"}),
        ]),

        html.Div([
            html.H4(id='button-clicks')]),
                dash_table.DataTable(id='DE-table',
                                     columns=[
                                         {'name': "Ensembl", "id": 'Ensembl'},
                                         {'name': "hgnc", "id": 'hgnc'},
                                         {'name': "wikigene_desc", "id": 'wikigene_desc'},
                                         {'name': "baseMean", "id": 'baseMean'},
                                         {'name': "log2FoldChange", "id": 'log2FoldChange'},
                                         {'name': "pvalue", "id": 'pvalue'},
                                         {'name': "padj", "id": 'padj'}],
                                     style_table={
                                         'maxHeight': '600px',
                                         'overflowY': 'scroll'
                                     },
                ),


        #dash_table.DataTable(
        #    id='DE-table',
        #    columns=[
        #            {'name': "Ensembl", "id": 'Ensembl'},
        #             {'name': "hgnc", "id": 'hgnc'},
        #             {'name': "wikigene_desc", "id": 'wikigene_desc'},
        #             {'name': "baseMean", "id": 'baseMean'},
        #             {'name': "log2FoldChange", "id": 'log2FoldChange'},
        #             {'name': "pvalue", "id": 'pvalue'},
        #             {'name': "padj", "id": 'padj'}],
        #    filter_action='custom',
        #    filter_query=''
        #),

        html.Br(),
        html.Div(id = 'table-box'),
        html.Div(dt.DataTable(id = 'de_table_comp', data=[{}]), style={'display': 'none'}),
        html.Div(id='de_table_comparison', style={'display': 'none'}),
        html.Br(),


        html.Br(),

        html.Div([
            dbc.Button('Enrichr', id='btn-enrichr', n_clicks=None),
        ]),

        html.Div([
            dbc.Tabs(
                [
                    dbc.Tab(label="GO: Biological process 2018 Up ", tab_id="tab-1"),
                    dbc.Tab(label="GO: Biological process Down", tab_id="tab-2"),
                    dbc.Tab(label="EGO: Cellular Component 2018 Up", tab_id="tab-3"),
                    dbc.Tab(label="GO: Cellular Component 2018 Down", tab_id="tab-4"),
                    dbc.Tab(label="GO: Molecular function 2018 Up", tab_id="tab-5"),
                    dbc.Tab(label="GO: Molecular function 2018 Down", tab_id="tab-6"),
                    dbc.Tab(label="Kegg 2016 Up", tab_id="tab-7"),
                    dbc.Tab(label="Kegg 2016 Down", tab_id="tab-8"),
                ]   ,
                id="enrichr_tabs",
                active_tab="tab-1",
            ),

            html.Div(id="enrichr_content"),
        ]),


        html.Div(id='Enrichr_GO_bp_up_ph', style={'display': 'none'}),
        html.Div(id='Enrichr_GO_bp_dn_ph', style={'display': 'none'}),
        html.Div(id='Enrichr_GO_cell_up_ph', style={'display': 'none'}),
        html.Div(id='Enrichr_GO_cell_dn_ph', style={'display': 'none'}),
        html.Div(id='Enrichr_GO_mf_up_ph', style={'display': 'none'}),
        html.Div(id='Enrichr_GO_mf_dn_ph', style={'display': 'none'}),
        html.Div(id='Enrichr_kegg_up_ph', style={'display': 'none'}),
        html.Div(id='Enrichr_kegg_dn_ph', style={'display': 'none'}),


    ], style={'padding': '25px'}

)
