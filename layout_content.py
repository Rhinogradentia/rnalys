
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
import base64

from dash import dcc
from dash import html
from dash import Input, Output, State
from dash import dash_table
import dash_daq as daq
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


def b64_image(image_filename):
    with open(image_filename, 'rb') as f:
        image = f.read()
    return 'data:image/png;base64,' + base64.b64encode(image).decode('utf-8')

#df_counts_combined_file = '/home/cfrisk/Dropbox/dash/data/df_counts_combined.tab'
#df_counts_combined = pd.read_csv(df_counts_combined_file, sep='\t', index_col=0)

#df_meta_combined_file = '/home/cfrisk/Dropbox/dash/data/df_meta_v3.tab'
#df_meta_combined = pd.read_csv(df_meta_combined_file, sep='\t', index_col=0)
#available_tissues = df_meta_combined['tissue'].unique()

start_genes = ['LEP', 'ADIPOQ', 'IL6', 'CCL2', 'RBP4', 'NAMPT', 'ITLN1']

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

dataset_file_path = 'data/datasets/datasets.csv'



layout_page1 = html.Div(style={'backgroundColor': '#f5f5f5','padding': '25px'},
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
                html.P('Load dataset:', style={'width': '100px', 'display': 'flex', 'verticalAlign': "middle", 'padding': '2px', 'margin-top': '15px'}),
                html.Div(id='dataset_loader_start', children=[
                    #html.Div(id='dataset_loader_start', style={'display': 'none'}),
                    dcc.Dropdown(id='datasets', value='New', placeholder='New'),
                ], style={'width': '400px', 'display': 'inline-block'}),
                html.Div([
                    html.I(className="fas fa-question-circle fa-lg", id="target", style={'padding': '10px'}),
                    dbc.Tooltip("Load dataset if you have any saved", target="target", style={'font-size': 15,'width': '400px', 'height': '20px', 'padding': '10px'}),
                ], style={'display': 'flex', 'align-items': 'center'}),
            ], style={'display': 'flex', 'align-items': 'center', 'width': '600px'}),
        ], style={'display': 'flex', 'align-items': 'center', 'width': '100%'}),

        html.P(id='dataset_name_placeholder'),

        #SELECT SAMPLES DIV
        html.Div([
            html.Div([
                html.Div([
                    html.H3('Select samples', style={'textAlign': 'left', 'padding': '5px', 'margin-top': '10px'}),

                    html.Div([
                        html.P('Variable 1:', style={'width': '120px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),

                        html.Div([
                            dcc.Dropdown(
                                id='variable1_selected_dropdown',
                                #options=[{'label': j, 'value': j} for j in df_meta_combined['type'].unique().tolist() + ['HF']],
                                #value=[''],
                                multi=True,
                            )
                        ], style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),

                    html.Div([
                        html.P('Variable 2:',style={'width': '120px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),
                        html.Div([
                            dcc.Dropdown(
                            id='variable2_selected_dropdown',
                            #value=[''],
                            multi=True)
                        ], style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),
                    #
                    html.Div([
                        html.P('Batch:', style={'width': '120px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),
                        html.Div([
                            dcc.Dropdown(
                            id='variable3_selected_dropdown',
                            #options=[{'label': j, 'value': j} for j in df_meta_combined['SeqTag'].unique()],
                            multi=True)
                        ], style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),

                    html.Div([
                        html.Div(id='exclude_list', style={'display': 'none'}),
                        html.P('Exclude:', style={'width': '120px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),
                        html.Div([
                            dcc.Dropdown(
                                id='exclude_dropdown',
                                #options=[{'label': j, 'value': j} for j in df_counts_combined.columns],
                                multi=True,
                                placeholder='Select Samples to exclude',

                            )
                        ],style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),
                    #
                    html.Div(id='save_dataset_div', children=[
                        dbc.Input(
                            id='save_dataset_name',
                            type='text',
                            placeholder='name dataset',
                            style={'padding': '5px', 'width': '300px'}
                        ),
                        html.Div([
                            dmc.Alert(
                                "-",
                                id="alert_datafiles_selection",
                                color='info',
                                withCloseButton=True,
                            ),
                        ], style={'display': 'none'}),  # style={'textAlign': 'left', 'width': '30%', 'margin-top': 5}),

                        dbc.Button('Save', id='btn_save_dataset', n_clicks=0, style={'margin-left': 10}),
                        html.Div(id='name_dataset_save'),
                    ], style={'width': '30%', 'display': 'flex', 'padding': '5px'}),

                    #html.Div(id='read_table'),
                    html.Div(id='selected_data', style={'display': 'none'}),
                    html.Div(id='intermediate-table', style={'display': 'none'}),

                ], className='row', style={'width': '40%','display': 'inline-block', 'margin-bottom': '10px'}),

                html.Div([
                    html.H3('Transformation & Normalization', style={'textAlign': 'left', 'padding': '10px', 'margin-top': '10px'}),

                    html.Div([
                           html.Div([
                               html.P('Select transf/norm:',
                                      style={'width': '300px', 'display': 'flex', 'verticalAlign': "middle",
                                             'padding': '2px'}),
                               html.Div([
                                   dcc.Dropdown(
                                       id='transformation',
                                       options=[{'label': j, 'value': j} for j in ['Sizefactor normalization', 'vsd','rlog','quantilenorm_log2']],
                                       multi=False,
                                       placeholder='Select transformation',)
                               ], style={'display': 'inline-block', 'width': '100%', 'verticalAlign': "middle"})
                           ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),
                            #], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),

                        html.Div([
                            html.P('Remove Confounding:', style={'width': '300px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'2px'}),
                            html.Div([
                                dcc.Dropdown(
                                    id='rm_confounding',
                                    #options=[{'label': j, 'value': j} for j in df_meta_combined.columns],
                                    multi=False,
                                    value=None)
                            ], style={'display': 'inline-block', 'width': '100%', 'verticalAlign': "middle"})
                        ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),



                        #)

                        html.Div([
                            daq.BooleanSwitch('force_run', on=False, color='#9B51E0', label=' Force re-run', labelPosition='left'),

                        ], style={'display':'inline-block', 'margin-top': '5px', 'margin-bottom': '5px'}),
                        html.Br(),

                        html.Div(
                            [
                                html.Div([
                                    dbc.Button('Submit', id='btn-selected_data_submit', n_clicks=0, size='md'),
                                ], style={'width':'30%' ,'display': 'flex', 'align-items': 'center', 'margin-top': '10px',
                                          'margin-right': '5px'}),

                                html.Div([
                                    dbc.Spinner(html.Div(id="submit-done"), size='md', delay_hide=1, delay_show=1,
                                                show_initially=False, color="primary"),
                                ], style={'display': 'flex', 'align-items': 'center', 'margin-left': '25px'}),

                            ],
                            style={'width': '80%', 'display': 'flex', 'flex-direction': 'row'}
                        ),
                    ]),
                ], className='row', style={'width': '50%', 'display': 'inline-block'})

            ], className='row', style={'top':'-25px', 'margin-left':'5px', 'margin-right':'5px', 'backgroundColor':'white','display': 'flex', 'border': '1px solid #C6CCD5', 'border-radius': '5px','justify-content': 'space-between'}),

            html.Div([
                dmc.Alert(
                    "-",
                    id="alert_selection",
                    color='info',
                    #is_open=True,
                    withCloseButton=True,
                    # n_clicks=0,
                ),
            ], style={'display':'none'}),#style={'textAlign': 'left', 'width': '30%', 'margin-top': 5}),

            html.Div([
                dmc.Alert(
                    "-",
                    id="alert_main_table",
                    color='info',
                    withCloseButton=True
                    #n_clicks=0,
                ),
            ], style={'display':'none'}),
        ], className='row', style={'z-index': '2', 'position': 'relative', 'display': 'flex', 'justify-content': 'space-between', 'margin-top': 5, 'margin-left': 5, 'margin-right': 5}),

        html.Div(id='select_to_explore_block',
                 style={
                     'width': '15px',  # width of the block
                     'height': '40px',  # height of the block
                     'backgroundColor': '#f5f5f5',
                     'z-index': '1',
                     'left': '20%',
                     #'top': '-10px',
                     'position': 'relative',
                 }
        ),

        html.Div([
            html.Div([
                html.Div([
                    html.H3('Explore genes', style={'textAlign': 'left', 'margin-top': 5}),

                    html.Div(id='Select_gene_input'),
                    html.Div([
                        html.Div([html.P('Ensembl:')], style={'width': '80px', 'display': 'inline-block'}),
                        html.Div([dcc.Dropdown(id='gene_list',#'input-2',
                                     #options=[{'label': j, 'value': j} for j in df_counts_combined.index]
                                     multi=True,
                                     value=['ENSG00000232810'])],
                                     style={'width': '70%', 'display': 'inline-block', 'verticalAlign': "middle"})
                    ], style={'backgroundColor':'white', 'margin-top': '10px'}),

                    html.Div([
                        html.Div([html.P('hgnc:')], style={'width': '80px', 'display': 'inline-block'}),
                        html.Div([dcc.Dropdown(id='input-2_hgnc', options=[{'label': j, 'value': j} for j in hgnc_dropdown],
                                     multi=True, value=start_genes)],style={'backgroundColor':'white','width': '70%', 'display': 'inline-block', 'verticalAlign':"middle"})
                    ]),
                ], style={'backgroundColor':'white','display': 'inline-block', 'width':'70%', 'verticalAlign':"middle"}),

                html.Div([

                    html.Div([
                        html.H4('Display options'),

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
                                          {'label': ' Full text', 'value': 'fulltext'},
                                      ], labelClassName='Full text'),

                    ], style={'backgroundColor':'white', 'verticalAlign':"middle"}

                    ),

                ], style={'margin-top': 10, 'border':
                    '1px solid #C6CCD5', 'padding': 15,
                          'border-radius': '5px', 'display':'inline-block', 'width': '30%', 'verticalAlign':"middle"}
                ),
            ], className='row', style={'width': '100%','display': 'inline-block', 'verticalAlign':"middle", 'margin-bottom': '10px'}),

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
            ], style={'margin-top': 5, 'backgroundColor':'white',}),

            #dcc.Graph(id='indicator-graphic2'),
            #dcc.Graph(id='hpa_graph'),
            #dash_table.DataTable(id='gene_table',
            #                     columns=[
            #                         {"name": i, "id": i} for i in sorted(df_symbol.columns)
            #                     ]),

        ],  className='row', style={'z-index': '2','position':'relative', 'backgroundColor':'white', 'display': 'flex', 'margin-left': '5px', 'margin-right':'5px', 'border': '1px solid #C6CCD5', 'padding': '20px', 'border-radius': '5px','justify-content': 'space-between'}),

        html.Div(id='sample_to_pca_block',
            style={
                'width': '15px',  # width of the block
                'height': '40px',  # height of the block
                'backgroundColor':  '#f5f5f5',
                'z-index': '1',
                'left': '20%',
                #'top': '-1px',
                'position': 'relative',

            }
        ),

        html.Div([
            html.Div([html.H2('PCA analysis')], style={"textAlign": "left", 'margin-top': 10, 'margin-left': 10}),


            #dbc.Input(id='input-gene', type='text', placeholder='Insert Gene (Ensembl)'),

            html.Div([
                dbc.Input(id='number_of_genes', type='text', placeholder='All Genes',
                          style={'width': '300px', 'height': '35px', 'padding': '20px', 'margin-right': '100px'}),
                html.P('Color by:',
                       style={'width': '100px', 'margin-left': '80px', 'margin-right': '10px', 'align-self': 'center',
                              'margin': '0', 'line-height': '35px'}),
                dcc.Dropdown(
                    id='meta_dropdown',
                    # options=[{'label': j, 'value': j} for j in df_meta_combined.columns],
                    value='tissue',
                    style={'width': '400px', 'textAlign': 'left', 'margin-right': '100px'}
                ),

                html.P('Group by:',
                       style={'width': '100px', 'margin-left': '100px', 'margin-right': '10px', 'align-self': 'center',
                              'margin': '0', 'line-height': '35px'}),
                dcc.Dropdown(
                    id='meta_dropdown_groupby',
                    # options=[{'label': j, 'value': j} for j in df_meta_combined.columns],
                    value='tissue',
                    style={'width': '400px', 'textAlign': 'left'}
                ),
                html.Div([
                    daq.BooleanSwitch(id='sample_names_toggle', on=False, label="Sample names",
                                      labelPosition="left"),
                ], style={'display': 'inline-block', 'margin-left': '100px'}),

            ], style={'display': 'flex', 'align-items': 'center', 'width': '100%', 'margin-top': 5}),

            dcc.Checklist(
                id='biplot_radio',
                options=[{'label': 'Biplot', 'value': 'biplot'}],
                style={'display':'none'}
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
                html.H4('PC correlation to variables'),

                html.Div([
                    html.I(className="fas fa-question-circle fa-lg", id="qMark_pc_cor", style={'padding': '10px', 'margin-bottom': '10px'}),
                    dbc.Tooltip("Columns with NaN values filtered", target="qMark_pc_cor", style={'font-size': 15,'width': '400px', 'height': '20px', 'padding': '10px'}),
                ], style={'display': 'flex', 'align-items': 'center'}),

            ], style={'display': 'flex', 'align-items': 'center', 'width': '600px'}),
            html.Div([
                #html.Div([dcc.Graph(id='clustergram')]),
                html.Div([dcc.Graph(id='pca_correlation')]),
            ]),

        ], className='row', style={'z-index': '2','position':'relative', 'backgroundColor':'white', 'display': 'flex', 'margin-left':'5px', 'border': '1px solid #C6CCD5', 'border-radius': '5px', 'margin-right':'5px'}),

        html.Div(id='pca_to_de_block',
                 style={
                     'width': '15px',  # width of the block
                     'height': '40px',  # height of the block
                     'backgroundColor': '#f5f5f5',
                     'z-index': '1',
                     'left': '20%',
                     # 'top': '-1px',
                     'position': 'relative',

                 }
        ),

        html.Div([
            html.Div([
                html.H3('DE analysis'),
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
                        html.P('Row sum count >:',
                               style={'width': '10%', 'display': 'flex', 'verticalAlign': "middle", 'padding': '5px'}),
                        html.Div([
                            dbc.Input(id='rowsum', type='text', value=10,
                                      placeholder='Select rowsum')
                        ], style={'display': 'inline-block', 'width': '20%', 'verticalAlign': "middle"})
                    ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),

                    html.Div([
                        html.P('Design formula:',
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
                        dbc.Button('Run analysis', id='btn-DE', n_clicks=None, style={'margin-right': '5px'}),
                        #], style={'width': '40%', 'padding': '5px'}),

                        dbc.Button('Export table', id='btn_export', style={'margin-right': '5px'}),
                        dbc.Button('Export plot', id='btn_export_plot', style={'margin-right': '5px'}),
                        html.Div(
                            dbc.Spinner(html.Div(id="submit-de-done"), size='md', delay_hide=1, delay_show=1,
                                        show_initially=False, color="primary"),
                            style={'align-self': 'center', 'margin-left': '25px'}
                        ),
                    ], style={'width': '100%', 'padding': '5px', 'display': 'flex', 'flex-direction': 'row', 'align-items': 'center'}
                ),


                #], style={'width': '100%', 'padding':'5px'})
            ], style={'margin-top': 10, 'padding': 5, 'backgroundColor':'white',}),


            html.Br(),



            html.P(id='export_placeholder'),
            html.Div(id='intermediate-DEtable', style={'display': 'none'}),
            html.Div(id='temp', style={'display': 'none'}),
            html.Div(id='pvalue', style={'display': 'none'}),


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
                ], style={'display': 'inline-block', 'verticalAlign':"middle", 'margin-left': '5px'}),


                html.Div([dbc.Button('Add DE set', id='add_de_table')],
                         style={'display': 'inline-block', 'verticalAlign':"middle", 'margin-left': '5px'}),
                html.Div([
                    dbc.Button('clear DE set', id='clear_de_table')],
                        style={'display': 'inline-block', 'verticalAlign':"middle", 'margin-left': '5px'}),
            ]),

            #TODO remake table to interactive table
            #TODO add Diverging Data Bars for log2Foldchange, maybe padj
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

            html.Br(),
            html.Div(id = 'table-box'),
            html.Div(dt.DataTable(id = 'de_table_comp', data=[{}]), style={'display': 'none'}),
            html.Div(id='de_table_comparison', style={'display': 'none'}),
            html.Br(),

        ], className='row', style={'backgroundColor':'white', 'display': 'flex', 'margin-left': '5px', 'margin-right':'5px', 'border': '1px solid #C6CCD5', 'padding': '20px', 'border-radius': '5px','justify-content': 'space-between'}),
        html.Br(),

        html.Div(id='de_to_enrichr_block',
                 style={
                     'width': '15px',  # width of the block
                     'height': '40px',  # height of the block
                     'backgroundColor': '#f5f5f5',
                     'z-index': '1',
                     'left': '20%',
                     # 'top': '-1px',
                     'position': 'relative',

                 }
        ),

        html.Div([

            html.Div(
                [
                    html.Div([
                        dbc.Button('Enrichr', id='btn-enrichr', n_clicks=None),
                    ], style={'width': '20%', 'display': 'flex', 'flex-direction': 'row', 'margin-top': '10px'}),
                    html.Div([
                        dbc.Spinner(html.Div(id="submit_done_enrichr"), size='md', delay_hide=1, delay_show=1,
                                    show_initially=False, color="primary"),
                    ], style={'margin-left': '25px', 'display': 'flex', 'align-items': 'center'}),
                ], style={'width': '40%', 'display': 'flex', 'flex-direction': 'row', 'align-items': 'center'}
            ),

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
            ], style={'margin-top': 5, 'padding': 5, 'display': 'inline-block'}),


            html.Div(id='Enrichr_GO_bp_up_ph', style={'display': 'none'}),
            html.Div(id='Enrichr_GO_bp_dn_ph', style={'display': 'none'}),
            html.Div(id='Enrichr_GO_cell_up_ph', style={'display': 'none'}),
            html.Div(id='Enrichr_GO_cell_dn_ph', style={'display': 'none'}),
            html.Div(id='Enrichr_GO_mf_up_ph', style={'display': 'none'}),
            html.Div(id='Enrichr_GO_mf_dn_ph', style={'display': 'none'}),
            html.Div(id='Enrichr_kegg_up_ph', style={'display': 'none'}),
            html.Div(id='Enrichr_kegg_dn_ph', style={'display': 'none'}),

        ], style={'backgroundColor':'white', 'display': 'flex', 'flexWrap': 'wrap', 'alignItems': 'flex-start', 'margin-left':'5px', 'margin-right':'5px', 'border': '1px solid #C6CCD5', 'padding': '20px', 'border-radius': '5px','justify-content': 'space-between'}),
    ]#, style={}

)

layout_index = html.Div([
        dbc.Container([
        
       
           

            html.Div([
                html.Div([
                    html.H5('Upload your counts with corresponding info file and get started!', style={'flex-grow': 1, 'margin-left': '10px'}),
                    html.Div([
                        html.I(className="fas fa-question-circle fa-lg", id="target", style={'padding': '10px'}),
                        dbc.Tooltip("Inserted count data should have columns matching index rows in meta data, the first column in the meta table should represent the sample", 
                                    target="target",
                                    style={'font-size': '15px'})
                    ], style={'display': 'inline-block', 'vertical-align': 'middle'}),
                ], style={'width': '70%', 'display': 'flex', 'align-items': 'center', 'justify-content': 'space-between'}),
            ]),


                html.Div([
                    html.Div([
                        html.Div(
                            children=[

                            dcc.Upload(
                                id='upload-data',
                                children=html.Div([
                                    'Drag and Drop or ',
                                    html.A('Select Count File')
                                ]),
                                style={
                                    'width': '400px',
                                    'height': '60px',
                                    'lineHeight': '60px',
                                    'borderWidth': '1px',
                                    'borderStyle': 'dashed',
                                    'borderRadius': '5px',
                                    'textAlign': 'center',
                                    'margin': '10px'
                                },
                                multiple=False
                            ),

                            html.Div(id='checkmark_counts_div',
                                children = [
                                html.Img(src=b64_image('assets/checkmark.jpg'), height='60',width='60'),
                            ], style={'display':'none'})

            ],style={'width':'100%', 'display': 'flex', 'justify-content': 'space-between', 'align-items': 'center'}),

                html.Div([
                    dcc.Upload(
                        id='upload-data-info',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Info File')
                        ]),
                        style={
                            'width': '400px',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '10px'
                        },
                        multiple=False
                    ),
                    html.Div(id='checkmark_info_div',
                             children=[
                                 html.Img(src=b64_image('assets/checkmark.jpg'), height='60',width='60'),
                             ], style={'display': 'none'})

                ], style={'width':'100%', 'display': 'flex', 'justify-content': 'space-between'}),
            ], style={'display': 'flex', 'justify-content': 'space-between', 'align-items': 'center'}),

        html.Div(id='alert_import_div', children=[
                dmc.Alert(
                    "import error",
                    id="alert_import",
                    color='info',
                    #is_open=True,
                    withCloseButton=True,
                    # n_clicks=0,
                ),
        ], style={'display':'none'}),

        html.Div(id='alert_import_info_div', children=[
            dmc.Alert(
                "Columns in Count data does not match Rows in Info data",
                id="alert_import_info",
                title='Import Error!',
                color='red',

                # is_open=True,
                withCloseButton=True,
                # n_clicks=0,
            ),
        ], style={'display':'none'}),

        html.Br(),

        html.Div(id='df_combability_check', style={'display':'none'}),
        #html.Div(id='df_info_check', style={'display':'none'}),

        html.Div(id='varaible_selection_div', children=[

            html.H5('Select variables for sample selection'),
            html.Div([
                html.P('Select variable 1 column', style={'width': '250px', 'display': 'flex', 'verticalAlign': "middle", 'padding': '2px'}),
                dcc.Dropdown(id='variable_selection1', style={'width': '250px', 'display': 'inline-block'}),
            ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%'}),

            html.Div([
                html.P('Select variable 2 column',
                       style={'width': '250px', 'display': 'flex', 'verticalAlign': "middle", 'padding': '2px'}),
                dcc.Dropdown(id='variable_selection2', style={'width': '250px', 'display': 'inline-block'}),
            ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%'}),

            html.Div([
                html.P('Select batch column',
                       style={'width': '250px', 'display': 'flex', 'verticalAlign': "middle", 'padding': '2px'}),
                dcc.Dropdown(id='variable_selection3', style={'width': '250px', 'display': 'inline-block'}),
            ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%'}),


        ], style={'display':'none'}),
        html.Br(),
        html.Div(id='page_proceed', children=[
            dbc.Button('Proceed', href='/page-1'),
        ], style={'display':'none'}),
        html.Br(),
        #dcc.Link('WGCNA + Enrichr', href='/page-2'),
        #html.Div(id='df_counts', style={'display': 'none'}),
        #html.Div(id='df_info', style={'display': 'none'}),
    ], style={'align-items': 'center', 'verticalAlign': 'center', 'align-content': 'center', 'justify-content': 'center', 'width': '100%', 'display': 'flex', 'flex-direction': 'column'}),
    ], style={'verticalAlign': 'center', 'border': '1px solid #C6CCD5', 'border-radius': '5px', 'padding': '50px', 'margin': '10%', }),
])
