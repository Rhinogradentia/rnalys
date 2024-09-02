import os
import json
import math
import flask
import dash
import dash_auth
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
from dash import dash_table as dt
import dash_daq as daq
from collections import Counter
from dash.exceptions import PreventUpdate
from sklearn.decomposition import PCA
from demos import dash_reusable_components as drc


#Genesymbol
df_symbol_file =  os.path.join('.', 'data', 'ensembl_symbol.csv')
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
df_volcano = pd.DataFrame(data, columns=columns)

datasets_txt_path = os.path.join('data', 'datasets.txt')

if os.path.isfile(datasets_txt_path):
    lPapers = []
    # Use the corrected path to open the file
    with open(datasets_txt_path, 'r') as f:
        for line in f:
            sPaper = line.split()[0]
            lPapers.append(sPaper)
    lPapers.append('New')
else:
    open(datasets_txt_path, 'a').close()


dataset_file_path = os.path.join('data', 'datasets', 'datasets.csv')

tabs_styles = {
    'height': '20px'
}

layout_page1 = html.Div(
        style={
        'backgroundColor': '#f5f5f5',
        'padding': '25px',
        'width': '100%',
        'margin': 'auto',
        'transform': 'scale(0.75)',  # Adjust this value to control the zoom level
        'transform-origin': 'top center'  # Adjusts the origin of transformation
    },
    children=[
        # Error Message
        html.Div(id="error-message"),
        # Top Banner
        html.Div(
            className="study-browser-banner row",
            children=[
                html.H2(className="h2-title", children="Rnalys"),
            ],
        ),

        html.Div([
            html.Div([
                html.P('Load dataset:', style={'width': '100px', 'display': 'flex', 'verticalAlign': "middle", 'padding': '2px', 'margin-top': '15px', 'display':'flex', 'flex-direction':'row'}),
                html.Div(id='dataset_loader_start', children=[
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
                    html.H5('Select samples', style={'textAlign': 'left', 'padding': '5px', 'margin-top': '10px'}),

                    html.Div([
                        html.P('Variable 1:', style={'width': '140px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),

                        html.Div([
                            dcc.Dropdown(
                                id='variable1_selected_dropdown',
                                multi=True,
                            )
                        ], style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),

                    html.Div([
                        html.P('Variable 2:',style={'width': '140px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),
                        html.Div([
                            dcc.Dropdown(
                            id='variable2_selected_dropdown',
                            multi=True)
                        ], style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),
                    #
                    html.Div([
                        html.P('Batch:', style={'width': '140px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),
                        html.Div([
                            dcc.Dropdown(
                            id='variable3_selected_dropdown',
                            multi=True)
                        ], style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),

                    html.Div([
                        html.Div(id='exclude_list', style={'display': 'none'}),
                        html.P('Exclude:', style={'width': '140px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'5px'}),
                        html.Div([
                            dcc.Dropdown(
                                id='exclude_dropdown',
                                multi=True,
                                placeholder='Select Samples to exclude',

                            )
                        ],style={'display': 'inline-block', 'width':'100%', 'verticalAlign':"middle"})
                    ], style={'display': 'flex', 'verticalAlign':"middle", 'width': '80%'}),
                    
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

                    html.Div(id='selected_data', style={'display': 'none'}),
                ], className='row', style={'width': '50%','display': 'inline-block', 'margin-bottom': '10px'}),

                html.Div([
                    html.H5('Transformation & Normalization', style={'textAlign': 'left', 'padding': '10px', 'margin-top': '10px'}),

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
                                       placeholder='Select transformation',
                                       value='vsd')
                               ], style={'display': 'inline-block', 'width': '100%', 'verticalAlign': "middle"})
                           ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),

                        html.Div([html.H6('Optional')

                        ]),

                        html.Div([
                            html.P('Remove Confounding:', style={'width': '300px', 'display': 'flex', 'verticalAlign':"middle", 'padding':'2px'}),
                            html.Div([
                                dcc.Dropdown(
                                    id='rm_confounding',
                                    #options=[{'label': j, 'value': j} for j in df_meta_combined.columns],
                                    multi=False,
                                    value=None)
                            ], style={'display': 'inline-block', 'width': '100%', 'verticalAlign': "midd    le"})
                        ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '80%'}),
 
                        html.Div(
                        [
                            # LEFT SIDE
                            html.Div(
                                    dbc.Switch(
                                        id="force_run",
                                        label="Force Re-run",
                                        value=False,
                                    ),
                                    style={
                                        'display': 'flex', 
                                        'align-items': 'left', 
                                        'justify-content': 'left',  
                                        'margin-top': '50px', 
                                        'margin-bottom': '5px',
                                        'width': '50%',  
                                    }
                            ),
                        
                            # RIGHT SIDE
                            html.Div(
                                [
                                    html.Div(

                                        children=[
                                            # Main Div container
                                            html.Div(
                                                children=[
                                            # Submit button
                                            dbc.Button(
                                                'Submit', 
                                                id='btn-selected_data_submit', 
                                                n_clicks=0, 
                                                size='md',
                                                style={'display': 'inline-block'}  # Add margin to the right of the button
                                            ),
                                            
                                            # Loading button with spinner (initially hidden)
                                            dcc.Loading(
                                                id="loading-indicator",
                                                type="circle", 
                                                children=[
                                                    html.Div(id='intermediate-table', style={'display': 'none'}),
                                                    html.Div(id="submit-done"),
                                                ],
                                                fullscreen=False,
                                                style={'margin-left':'50px'}
                                            )
                                        ],
                                        style={
                                            'display': 'flex',
                                            'align-items': 'center',
                                            'justify-content': 'flex-start', 
                                            'width': '100%',
                                            'margin-top':'50px'
                                        }
                                            ),
                                        ]
                                )
                                ], 
                                style={
                                    'display': 'flex', 
                                    'align-items': 'center',  # Center align vertically
                                    'justify-content': 'center',  # Center horizontally if needed
                                    'width': '50%'  # Adjust width as necessary
                                }
                            ),
                        ],
                        style={
                            'display': 'flex', 
                            'align-items': 'center',  # This ensures both sides are at the same level
                            'width': '100%',  # Ensure parent fills out the width
                        }
                    )
                    ]),
                ], className='row', style={'width': '50%', 'display': 'inline-block'})

            ], className='row', style={'top':'-25px', 'margin-left':'5px', 'margin-right':'5px', 'backgroundColor':'white','display': 'flex', 'border': '1px solid #C6CCD5', 'border-radius': '5px','justify-content': 'space-between'}),

            html.Div([
                dmc.Alert(
                    "-",
                    id="alert_selection",
                    color='Warning',
                    #is_open=True,
                    withCloseButton=True,
                    # n_clicks=0,
                ),
            ], style={'display':'none'}),#style={'textAlign': 'left', 'width': '30%', 'margin-top': 5}),

            html.Div([
                dmc.Alert(
                    "Select transformation",
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
                    html.H5('Explore genes', style={'textAlign': 'left'}),

                    html.Div(id='Select_gene_input'),
                    html.Div([
                        html.Div([html.P('Ensembl:')], style={'width': '80px', 'display': 'inline-block'}),
                        html.Div([dcc.Dropdown(
                                    id='gene_list',
                                    multi=True,
                                     value=[])],
                                     style={'width': '70%', 'display': 'inline-block', 'verticalAlign': "middle"})
                    ], style={'backgroundColor':'white', 'margin-top': '10px'}),

                    html.Div([
                        html.Div([html.P('hgnc:')], style={'width': '80px', 'display': 'inline-block'}),
                        html.Div([dcc.Dropdown(id='input-2_hgnc', options=[{'label': j, 'value': j} for j in hgnc_dropdown],
                                     multi=True, value=[])],style={'backgroundColor':'white','width': '70%', 'display': 'inline-block', 'verticalAlign':"middle"})
                    ]),
                ], style={'backgroundColor':'white','display': 'inline-block', 'width':'70%', 'verticalAlign':"middle"}),

                html.Div([

                    html.Div([
                        html.H6('Display options', style={'padding': '0px', 'margin-top': '1px', 'margin-bottom':'2px'}),

                    
                        html.Div([
                            html.P('Gene ID: ',
                                style={'width': '100px', 'display': 'flex', 'verticalAlign': "middle",
                                        'padding': '2px', 'margin-top':'3px'}),
                            html.Div([dcc.Dropdown(id='radio_symbol',
                                        options=[{'label': j+' ', 'value': j} for j in ['Ensembl ', 'Symbol']],
                                        multi=False,
                                        value=['Ensembl'])
                                    ], style={'display': 'inline-block', 'width': '80%', 'verticalAlign': "middle", 'margin-top':'5px'}),
                        ],  style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%'}),

                        html.Div([
                            html.P('Group by: ',
                                style={'width': '100px', 'display': 'flex', 'verticalAlign': "middle",
                                        'padding': '2px'}),
                            html.Div([dcc.Dropdown(id='radio-grouping',
                                            options=[
                                               {'label': 'Variable 1 ', 'value': 'var1'},
                                               {'label': 'Variable 2 ', 'value': 'var2'},
                                               {'label': 'Batch ', 'value': 'var3'},
                                           ],
                                        multi=False,
                                        value='var1',
                                        style={'height':'2px'})
                                    ], style={'display': 'inline-block', 'width': '80%', 'verticalAlign': "middle", 'height': '2px'}),
                        ],  style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%'}),

                        dcc.Checklist(id='full_text',
                                      options=[
                                          {'label': ' Full text', 'value': 'fulltext'},
                                      ], labelClassName='Full text'),

                    ], style={'backgroundColor':'white', 'verticalAlign':"middle"}

                    ),

                ], style={'margin-top': '1%', 'border':
                    '1px solid #C6CCD5', 'padding': 15,
                          'border-radius': '5px', 'display':'inline-block', 'width': '30%', 'verticalAlign':"middle"}
                ),
            ], className='row', style={'width': '100%','display': 'inline-block', 'verticalAlign':"middle", 'margin-bottom': '5px'}),

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
                        style={'height':'5vh'}
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
            html.Div([html.H5('PC Analysis')], style={"textAlign": "left", 'margin-top': 10, 'margin-left': 10}),
            html.Div([
                dbc.Input(
                    id='number_of_genes', type='text', placeholder='All Genes',
                    style={'height': '35px', 'width': '20%', 'margin-right': '2%'}
                ),
                html.Div(
                    children=[
                        html.P(
                            'Color by:',
                            style={'margin-right': '1%', 'line-height': '35px', 'margin-bottom': '0', 'margin-top': '4px'}
                        ),
                        dcc.Dropdown(
                            id='meta_dropdown',
                            style={'width': '70%'}
                        ),
                    ],
                    style={'display': 'flex', 'flex-direction': 'row', 'width': '25%', 'align-items': 'center'}
                ),
                html.Div(
                    children=[
                        html.P(
                            'Group by:',
                            style={'margin-right': '1%', 'line-height': '35px', 'margin-bottom': '0', 'margin-top': '4px'}
                        ),
                        dcc.Dropdown(
                            id='meta_dropdown_groupby',
                            style={'width': '70%'}
                        ),
                    ],
                    style={'display': 'flex', 'flex-direction': 'row', 'width': '25%', 'align-items': 'center'}
                ),
                daq.BooleanSwitch(
                    id='sample_names_toggle', on=False, label="Sample names", labelPosition="left",
                    style={'width': '10%', 'margin-left': '2%'}
                ),
            ], style={'display': 'flex', 'flex-direction': 'row', 'width': '100%', 'align-items': 'center'}),

            dcc.Checklist(
                id='biplot_radio',
                options=[{'label': 'Biplot', 'value': 'biplot'}],
                style={'display':'none'}
            ),

            dcc.RadioItems(id='biplot_text_radio'),

            html.Div(
    [
        # Parent div for side by side tabs layout
        html.Div(
            [
                # Left Side: Container for the PCA plots
                html.Div(
                    [
                        dbc.Tabs(
                            [
                                # Tab 1: Contains the 2D PCA plot
                                dbc.Tab(
                                    dcc.Graph(id='pca_and_barplot', style={'width': '100%'}),
                                    label="2 Components"
                                ),
                                # Tab 2: Contains the 3D PCA plot
                                dbc.Tab(
                                    dcc.Graph(id='pca-graph', style={'width': '100%'}),
                                    label="3 Components"
                                ),
                            ],
                            style={'width': '100%'}
                        )
                    ],
                    style={'width': '60%', 'display': 'inline-block'}
                ),

                # Right Side: Container for the Barplot and PCA Component Explained

                html.Div(
                    [
                        dbc.Tabs(
                            [
                                # Tab 1: Contains the barplot
                                dbc.Tab(
                                    dcc.Graph(id='barplot', style={'width': '100%'}),
                                    label="Bar Plot"
                                ),
                                # Tab 2: Contains the PCA Component Explained
                                dbc.Tab(
                                    dcc.Graph(id='pca_comp_explained', style={'width': '100%'}),
                                    label="PC Variance"
                                ),
                                # Tab 3:
                                dbc.Tab(
                                    dcc.Graph(id='pca_correlation', style={'width': '100%'}),
                                    label="PC correlations"
                                ),
                            
                            ],
                            style={'width': '100%'}
                        )
                    ],
                    style={'width': '40%', 'display': 'inline-block'}
                ),
            ],
            style={'display': 'flex', 'margin-top':'10px'}
            )
            ]
        ),

            html.Div(id='clicked'),
    

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
                html.H5('DE analysis'),
                 # Main container for side-by-side layout
                html.Div([
                # Left side container
                html.Div([
                    # Contains all left side input controls and buttons
                    html.Div([
                        # First section
                        html.Div([
                            html.Br(),
                            html.P('Program:',
                                style={'width': '140px', 'display': 'inline-block', 'verticalAlign': "middle", 'padding': '5px'}),

                            dcc.Dropdown(
                                id='program',
                                options=[{'label': j, 'value': j} for j in ['DESeq2', 'edgeR', 'limma']],
                                value='DESeq2',
                                multi=False,
                                style={'height': '35px', 'width': '70%', 'display': 'inline-block', 'verticalAlign': "middle"}
                            )
                        ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '10px', 'margin-top':'240px'}),

                        # Second section
                        html.Div([
                            html.P('Row sum >:',
                                style={'width': '140px', 'display': 'inline-block', 'verticalAlign': "middle", 'padding': '5px'}),
                            dbc.Input(id='rowsum', type='text', value=10,
                                    placeholder='Select rowsum',
                                    #sstyle={, 'padding': '5px'}
                                    style={'height': '35px','width': '70%', 'display': 'inline-block', 'verticalAlign': "middle"})
                        ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '10px'}),

                        # Third section
                        html.Div([
                            html.P('Design:',
                                style={'width': '140px', 'display': 'inline-block', 'verticalAlign': "middle", 'padding': '5px'}),
                            dbc.Input(id='design', type='text', value='~',
                                    placeholder='Design formula',
                                    style={'height': '35px', 'width': '70%', 'display': 'inline-block', 'verticalAlign': "middle"})
                        ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '10px'}),

                        # Fourth section
                        html.Div([
                            html.P('Reference:',
                                style={'width': '140px', 'display': 'inline-block', 'verticalAlign': "middle", 'padding': '5px'}),
                            dbc.Input(id='reference', type='text',
                                    placeholder='Enter reference for DE',
                                    style={'height': '35px','width': '70%', 'display': 'inline-block', 'verticalAlign': "middle"})
                        ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '10px'}),

                        # X, Y selection for plot
                        html.Div([
                            html.P('Display, X, Y',
                                style={'width': '140px', 'display': 'inline-block', 'verticalAlign': "middle", 'padding': '5px'}),

                            dcc.Dropdown(
                                id='volcano_xaxis',
                                options=[{'label': j, 'value': j} for j in ['log2FoldChange', 'baseMean']],
                                value='log2FoldChange',
                                multi=False,
                                placeholder='x-axis',
                                style={'height': '35px', 'width': '80%', 'display': 'inline-block', 'verticalAlign': "middle"}
                            ),

                            dcc.Dropdown(
                                id='volcano_yaxis',
                                options=[{'label': j, 'value': j} for j in ['log2FoldChange', '-log10(p)', 'baseMean']],
                                value='-log10(p)',
                                multi=False,
                                placeholder='y-axis',
                                style={'height': '35px', 'width': '80%', 'display': 'inline-block', 'verticalAlign': "left"}
                            )
                        ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '1px'}),
                    ], style={'width': '100%'}),

                    # Buttons
                    html.Div([
                        html.Div(
                            children=[
                                dbc.Button('Run analysis', id='btn-DE', n_clicks=None, style={'width': '49%', 'margin-right': '5px'}),
                                                
                                # Loading button with spinner (initially hidden)
                                dcc.Loading(
                                    id="loading-indicator",
                                    type="circle", 
                                    children=[
                                        html.Div(id='intermediate-DEtable', style={'display': 'none'}),
                                        html.Div(id="submit-de-done")
                                    ],
                                    fullscreen=False,
                                    style={'margin-left':'50px'},
                                )
                        ], style={'width': '100%', 'padding': '5px', 'display': 'flex', 'flex-direction': 'row', 'align-items': 'center'}),
                    ]),
                    html.Div([    
                        
                        dbc.Button('Export table', id='btn_export', style={'width': '80%', 'margin-right': '5px'}),
                        dbc.Button('Export plot', id='btn_export_plot', style={'width': '80%', 'margin-right': '5px'}),
                        
                    ], style={'width': '100%', 'padding': '5px', 'display': 'flex', 'flex-direction': 'row', 'align-items': 'center'})
                ], style={'width': '30%', 'display': 'inline-block', 'padding': '10px'}),

                    # Right side container


                    html.Div([
                        dbc.Tabs([
                            dbc.Tab(label='Volcano Plot', children=[
                                html.Div([
                                    html.Div(
                                        dcc.RangeSlider(
                                            id='volcanoplot-input',
                                            min=-7,
                                            max=7,
                                            step=0.05,
                                            marks={i: {'label': str(i)} for i in range(-8, 9)},
                                            value=[-1, 1]
                                        ),
                                        #style={'width': '100%', 'margin': 'auto'}#, 'padding': '20px'}
                                    ),
                                html.Div(
                                    dcc.Graph(
                                        id='volcanoplot',
                                        figure=dashbio.VolcanoPlot(
                                            effect_size='log2FoldChange', 
                                            p='padj', 
                                            gene='Ensembl',
                                            logp=True, 
                                            snp='Ensembl', 
                                            xlabel='log2FoldChange',
                                            genomewideline_value=2.5, 
                                            dataframe=df_volcano
                                        ),
                                        style={'height': '650px', 'width': '800px'}
                                    ),
                                )
                        ], style={'width': '100%', 'display': 'inline-block'})#style={'width': '100%', 'display': 'inline-block', 'verticalAlign': "top", 'padding': '10px'})
                        ],), #style={'height':'5vh', 'line-height': '5vh'}),

                            dbc.Tab(
                                        dcc.Graph(id='ma_plot', style={'width': '100%', 'margin': 'auto', 'height': '600px', 'width': '800px'}),
                                                    label="MA-plot"
                                                ),
                                    
                                ],)# style=tabs_styles)#style={'height': '20px'})
                        ])#, style={'width': '70%', 'display': 'inline-block', 'padding': '10px', 'height':'1%'})
                
                    ], style={'display': 'flex', 'width': '100%'}),

                    html.Br(),
                    html.P(id='export_placeholder'),
                    
                    html.Div(id='temp', style={'display': 'none'}),
                    html.Div(id='pvalue', style={'display': 'none'}),
                    html.Div(id='export_plot_clicked'),

        ], style={'margin-top': 30, 'padding': 5, 'backgroundColor': 'white'}),

        html.Div(id='number_of_degenes'),

        html.Div([
                
        html.Div([
        html.P('Filter on'),
        html.Div([
            
        #html.Br(),
        # Left Section: Filter on significance
        html.Div([
            
            html.P('Significance: ',
                   #style={'width': '200px', 'display': 'flex', 'verticalAlign': "middle", 'padding': '0px'}),
                   style={'width': '120px', 'display': 'inline-block', 'verticalAlign': "middle", 'padding': '5px'}),
            html.Div([
                dbc.Input(id='toggle_sig', type='text', value='0.05',
                          placeholder='Select p-value', 
                          style={'height': '35px', 'padding': '5px'})
            ],  style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '10px'})
        ],  style={'display': 'flex', 'verticalAlign': "middle", 'width': '20%', 'margin-bottom': '10px'}),#style={'display': 'flex', 'verticalAlign': "middle", 'alignItems': 'center', 'padding': '0px'}),

        # Spacer between sections
        html.Div(style={'width': '10px'}),  # Adjust or remove this spacer to change distance

        # Right Section: Basemean
        html.Div([
            html.P('Basemean: ',
                   style={'width': '100px', 'display': 'flex', 'verticalAlign': "middle", 'padding': '5px'}),
            html.Div([
                dbc.Input(id='toggle_basemean', type='text', value='0',
                          placeholder='basemean',
                          style={'height': '35px', 'padding': '5px'})
            ], style={'display': 'inline-block', 'width': '60%', 'verticalAlign': "middle"})
        ],  style={'display': 'flex', 'verticalAlign': "middle", 'width': '20%', 'margin-bottom': '10px'}),
        ], style={'display': 'flex', 'alignItems': 'center'}),

        ], style={'padding': '10px'}),

                html.Div([
                    dbc.Button(id='sig_submit', n_clicks=0, children='Submit'),
                ], style={'display': 'inline-block', 'verticalAlign':"middle", 'margin-left': '5px'}),


                html.Div([dbc.Button('Add DE set', id='add_de_table')],
                         style={'display': 'inline-block', 'verticalAlign':"middle", 'margin-left': '5px'}),
                html.Div([
                    dbc.Button('clear DE set', id='clear_de_table')],
                        style={'display': 'inline-block', 'verticalAlign':"middle", 'margin-left': '5px'}),
            ]),

            html.Div([
                html.H5(id='button-clicks')]),
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

            html.H5('Enrichr analysis', style={'textAlign': 'left', 'padding': '5px', 'margin-top': '10px'}),
            html.Div(
                [
                html.Div([    
                
                    dbc.Button('Run Enrichr', id='btn-enrichr', n_clicks=None, style={'width':'20%', 'height':'10%'}),
                    dbc.Button('Export plot', id='btn_export_enrichr_plot', n_clicks=None, style={'width': '20%', 'margin-left': '5px', 'height':'10%'}),
                    html.P(id='export_enrichr_plot'),
                    
                ], style={'width': '100%', 'display': 'flex', 'flex-direction': 'row'}),

                    
                ], style={'width': '60%', 'display': 'flex', 'flex-direction': 'row', 'align-items': 'center'}
            ),
            
            html.Div([
                dbc.Tabs(
                    [
                        dbc.Tab(label="GO: Biological process 2018 Up ", tab_id="tab-1"),
                        dbc.Tab(label="GO: Biological process Down", tab_id="tab-2"),
                        dbc.Tab(label="GO: Cellular Component 2018 Up", tab_id="tab-3"),
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

        ], className='row', style={'height':'1000','backgroundColor':'white', 'display': 'flex', 'flexWrap': 'wrap', 'alignItems': 'flex-start', 'margin-left':'5px', 'margin-right':'5px', 'border': '1px solid #C6CCD5', 'padding': '5px', 'border-radius': '5px','justify-content': 'space-between'}),
    ]#, style={}

)

layout_index = html.Div(
    style={
        'backgroundColor': '#f5f5f5',
        'padding': '25px',
        'width': '100%',
        'margin': 'auto',
        'transform': 'scale(0.8)',  # Adjust this value to control the zoom level
        'transform-origin': 'top center'  # Adjusts the origin of transformation
    },
    children=[
        dbc.Container([
            # Flex container for both left and right sides
            html.Div([
                # Left side content
                html.Div([
                    html.H5('Upload counts & info file to get started!', style={'margin-left': '10px', 'margin-bottom': '20px'}),

                    # First upload and its alert
                    html.Div([
                        # Container for the upload and checkmark
                        html.Div([
                            dcc.Upload(
                                id='upload-data',
                                children=html.Div([
                                    'Drag and Drop or select count file'
                                ], style={
                                    # Ensure the div inside takes full size of its parent to center text correctly
                                    'width': '100%',
                                    'height': '100%',
                                    'display': 'flex',
                                    'alignItems': 'center',  # Vertical center
                                    'justifyContent': 'center'  # Horizontal center
                                }),
                                style={
                                    'width': '300px',
                                    'height': '60px',
                                    'lineHeight': '60px',
                                    'borderWidth': '2px',
                                    'borderStyle': 'dashed',
                                    'borderRadius': '5px',
                                    'textAlign': 'center',
                                    'margin': '10px',
                                    'padding': '10px'
                                },
                                multiple=False
                            ),
                            html.Div(id='checkmark_counts_div',
                                children=[
                                    # Assuming you have a function `b64_image` to convert image path to base64 string
                                    html.Img(src=b64_image('assets/checkmark.jpg'), height='60', width='60'),
                                ], style={'display':'none'})
                        ], style={'display': 'flex',  'width': '100%'}),

                        # Alert next to the upload
                        html.Div(id='alert_import_div', children=[
                            dbc.Alert(
                                "Import error",
                                id="alert_import",
                                color='info',
                                dismissable=True,
                            ),
                        ], style={'display':'none', 'flex': '1'}),

                    ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '10px', 'margin-top':'20'}),

                    # Second upload and its alert
                    html.Div([
                        # Container for the upload and checkmark
                        html.Div([
                            dcc.Upload(
                                id='upload-data-info',
                                children=html.Div([
                                    'Drag and Drop or select info file'
                                ], style={
                                    # Ensure the div inside takes full size of its parent to center text correctly
                                    'width': '100%',
                                    'height': '100%',
                                    'display': 'flex',
                                    'alignItems': 'center',  # Vertical center
                                    'justifyContent': 'center'  # Horizontal center
                                }),
                                style={
                                    'width': '300px',
                                    'height': '60px',
                                    'lineHeight': '60px',
                                    'borderWidth': '2px',
                                    'borderStyle': 'dashed',
                                    'borderRadius': '5px',
                                    'textAlign': 'center',
                                    'margin': '10px',
                                    'padding': '10px'
                                },
                                multiple=False
                            ),
                            html.Div(id='checkmark_info_div',
                                     children=[
                                         html.Img(src=b64_image('assets/checkmark.jpg'), height='60', width='60'),
                                     ], style={'display': 'none'})

                        ], style={'display': 'flex',  'width': '100%'}),

                        html.Div(id='alert_import_info_div2', children=[
                            dbc.Alert(
                                "Import error",
                                id="alert_import_info_2",
                                color='info',
                                dismissable=True,
                            ),
                        ], style={'display':'none'}),

                    ], style={'display': 'flex', 'width': '100%', 'align-items': 'center'}),

                    html.Div([
                        html.H5('Or try:', style={'margin-left': '17px', 'margin-bottom': '20px', 'margin-top':'10px'}),

                    ]),

                    html.Div([
                        dbc.Button('Example data', id='generate-example-data', className='btn btn-dark'),
                        html.Div(id='example-data-output')  # Placeholder for displaying an alert or a link to the generated data
                    ], style={'text-align': 'left', 'width': '50%', 'display': 'inline-block', 'margin-left': '10px'}),

                ], style={'display': 'flex', 'flex-direction': 'column', 'width': '50%', 'padding-right': '10px', 'align-items': 'flex-start'}),

                # Right side content
                html.Div(id='variable_selection_div', children=[
                    html.H5('Select variables for sample selection', style={'margin-bottom': '20px'}),
                    html.Div([
                        html.P('Select variable 1 column', style={'width': '250px', 'verticalAlign': "middle", 'padding': '2px'}),
                        dcc.Dropdown(id='variable_selection1', style={'width': '250px'}),
                    ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '10px'}),

                    html.Div([
                        html.P('Select variable 2 column', style={'width': '250px', 'verticalAlign': "middle", 'padding': '2px'}),
                        dcc.Dropdown(id='variable_selection2', style={'width': '250px'}),
                    ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '10px'}),

                    html.Div([
                        html.P('Select batch column', style={'width': '250px', 'verticalAlign': "middle", 'padding': '2px'}),
                        dcc.Dropdown(id='variable_selection3', style={'width': '250px'}),
                    ], style={'display': 'flex', 'verticalAlign': "middle", 'width': '100%', 'margin-bottom': '10px'}),

                    html.Div([
                        html.Div(id='page_proceed', children=[
                            dbc.Button('Proceed', color='success', href='/page-1'),
                        ], style={'display':'none'})
                    ], style={'text-align': 'right', 'width': '93%', 'display': 'inline-block', 'margin-top':'50px'}),

                ], style={'display':'block', 'width': '50%', 'padding-left': '10px'}),

            ], style={'display': 'flex', 'width': '100%', 'align-items': 'flex-start'}),



            html.Div(id='alert_import_info_div', children=[
                dbc.Alert(
                    "Columns in Count data does not match Rows in Info data",
                    id="alert_import_info",
                    color='danger',
                    dismissable=True,
                ),
            ], style={'display':'none', 'flex': '1'}),

        ], style={'border': '1px solid #C6CCD5', 'border-radius': '5px', 'padding': '20px', 'margin': '2%', 'backgroundColor':'white'}),
    ]
)
