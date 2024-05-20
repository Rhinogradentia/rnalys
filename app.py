# -*- coding: utf-8 -*-

"""
Application Name: Rnalys
Description: This application is a Dash-based web app designed to quick and easy get analyse rna-seq data 
Author: Christoffer Frisk
Contact: Christofffer@gmail.com
Created on: 2022.01.01
Last Modified: []
License: MIT license

This application utilizes Flask and Dash to create a web-based interface for [describe the general use case]. 
It incorporates features such as [mention any significant features like data visualization, user interaction, etc.].

Modules and Libraries:
- Flask: Used for setting up the server.
- Dash: Main framework for the web application.
- Pandas, Numpy: Data manipulation and numerical calculations.
- Plotly: Creating interactive plots and graphs.
- [Other Libraries]: [Their purposes].
"""

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
import base64
import io
import random
import string
from typing import List, Tuple, Any
from io import StringIO
from dash import dcc
from dash import html
from dash import dash_table
from collections import Counter
from dash.exceptions import PreventUpdate
from sklearn.decomposition import PCA
from demos import dash_reusable_components as drc
from layout_content import layout_index, layout_page1
from dash.dependencies import Input, Output, State
from datetime import datetime

# Configuration for external stylesheets
FONT_AWESOME = "https://use.fontawesome.com/releases/v5.10.2/css/all.css"


app = dash.Dash(
    __name__,
    # external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css', dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
    # external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'],
    external_stylesheets=[dbc.themes.LITERA, FONT_AWESOME],
    suppress_callback_exceptions=True
)
"""
Initializes the Dash app with specific external stylesheets.
- dbc.themes.LITERA: Bootstrap theme for styling.
- FONT_AWESOME: Font Awesome for icons.
suppress_callback_exceptions is set to True for handling exceptions in callbacks.
"""

#df meta : df_meta_v3.tab
#df counts: df_counts_combined.tab

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),

    # Data shared between pages
    dcc.Store(id='df_counts'),  # , storage_type='session'),
    dcc.Store(id='df_info'),  # storage_type='session'),
    dcc.Store(id='variable_selection1_store'),
    dcc.Store(id='variable_selection2_store'),
    dcc.Store(id='variable_selection3_store'),

    # Spliting dataframes to speed up for referencing
    dcc.Store(id='counts_index_store'),
    dcc.Store(id='counts_columns_store'),
    dcc.Store(id='info_columns_store'),
    dcc.Store(id='info_index_store'),

    html.Div(id='page-content')
])

# Genesymbol
df_symbol_file = os.path.join('data', 'ensembl_symbol.csv')
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

lExclude = []
lTissue = []
PAGE_SIZE = 10

# Load Human protein atlas (consensus rna tissue)
df_hpa_file = os.path.join('data', 'rna_tissue_consensus.tsv')
df_hpa =  pd.read_csv(df_hpa_file, sep='\t', index_col='Gene')

# Setup for Enrichr
databases_enrichr = ['TRANSFAC_and_JASPAR_PWMs', 'OMIM_Disease', 'OMIM_Expanded', 'KEGG_2016', 'BioCarta_2016',
                     'WikiPathways_2016', 'Panther_2016', 'GO_Biological_process_2018', 'GO_Cellular_Component_2018',
                     'GO_Molecular_Function_2018', 'Human_Phenotype_Ontology', 'MGI_Mammalian_Phenotype_2017',
                     'Jensen_DISEASES']

# Create temporary df for volcanoplot init
columns = ['Ensembl', 'hgnc', 'baseMean', 'log2FoldChange', 'pvalue', 'padj']
data = [['ens', 'hgnc1', 10, 4, 0.005, 0.9995]]
df = pd.DataFrame(data, columns=columns)

# Layouts
layout_index = layout_index
layout_page1 = layout_page1

def generate_random_string(length=7):
    characters = string.ascii_letters + string.digits
    random_string = ''.join(random.choices(characters, k=length))
    return random_string


def write_dataset(lSamples, lExclude, id_name=None):

    """
    Writes dataset information to a CSV file. This function is used for storing sample names,
    excluded names, and an optional ID name.

    Parameters:
    - lSamples (list): List of sample names.
    - lExclude (list): List of names to exclude.
    - id_name (str, optional): An identifier name for the dataset.

    The function checks if the 'data/datasets/' directory exists, creates it if not, and then writes 
    the dataset information to 'data/datasets/datasets.csv'.
    """
    # Generate a random filename
    dataset_file_path = os.path.join('data', 'datasets', 'datasets.csv')
    sample_names = ' '.join(lSamples)
    exclude_names = ' '.join(lExclude)
    outdir = os.path.join('data', 'datasets')

    if not os.path.exists(outdir):
        # If the directory doesn't exist, create it
        os.makedirs(outdir)

    df_new = pd.DataFrame({
        'Sample Names': [sample_names],
        'Exclude Names': [exclude_names],
        'ID Name': [id_name]
    })

    # If the file exists, read the existing data and append the new data
    if os.path.isfile(dataset_file_path):
        df_existing = pd.read_csv(dataset_file_path)
        df = pd.concat([df_existing, df_new], ignore_index=True)
    else:
        df = df_new

    # Write the DataFrame to a CSV file
    df.to_csv(dataset_file_path, index=False)

def load_dataset(dataset_name):
    dataset_file_path = os.path.join('data', 'datasets', 'datasets.csv')

    if os.path.isfile(dataset_file_path):
        df = pd.read_csv(dataset_file_path, index_col=False)
        matching_row = df[df['ID Name'] == dataset_name]

        # If there is a matching row, return the Sample Names
        if not matching_row.empty:
            # print(matching_row['Sample Names'].values[0])
            return matching_row['Sample Names'].values[0], matching_row['Exclude Names'].values[0]

        # If there is no matching row, return None
        return None

    else:
        df = pd.DataFrame(columns=['Sample Names', 'Exclude Names', 'ID Name'])
        # Write the DataFrame to 'data/datasets/datasets.csv'
        df.to_csv(dataset_file_path, index=False)
        return None


def write_session_to_file(sample_names, rm_confounding, transformation, file_string, id_name=None):
    '''
    Writes session information to a file. This is used for tracking various parameters of a user session.
    
    Parameters:
        - sample_names (list): List of sample names.
        - rm_confounding (bool): Flag for removing confounding variables.
        - transformation (str): The type of transformation applied.
        - file_string (str): A string identifier for the file.
        - id_name (str, optional): An identifier name for the session.

        This function writes the session data to 'data/datasets/session_file.txt'.
    '''

    session_file_path = os.path.join('data', 'datasets', 'session_file.txt')
    # Write the sample names to the file
    # > sample1 sample2 sample3 transformation %id_name %file_string
    df = pd.DataFrame({
        'Sample Names': [' '.join(sample_names)],
        'Transformation': [transformation],
        'ID Name': [id_name],
        'rm_confounding': [rm_confounding],
        'file_string': [file_string]
    })

    if os.path.isfile(session_file_path):
        df_existing = pd.read_csv(session_file_path, encoding='utf-8')
        df = pd.concat([df_existing[df.columns.tolist()], df], ignore_index=True)
        df.to_csv(session_file_path, index=False, encoding='utf-8')
    else:
        df.to_csv(session_file_path, index=False, encoding='utf-8')

def search_session(sample_string, transformation, rm_confounding):
    """
    Searches for a session in the session file that matches given parameters.

    Parameters:
    - sample_string (list): List of sample names.
    - transformation (str): The type of transformation applied.
    - rm_confounding (bool): Flag for removing confounding variables.

    Returns:
    - list or None: Returns a list of file strings if a matching session is found, otherwise None.
    """

    session_file_path = os.path.join('data', 'datasets', 'session_file.txt')
    outdir = os.path.join('data', 'datasets')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if not os.path.isfile(session_file_path):
        df = pd.DataFrame(columns=['Sample Names', 'ID Name', 'Transformation', 'rm_confounding', 'file_string'])
        df.to_csv(session_file_path, index=False)

    df = pd.read_csv(session_file_path)

    # Check if the given values match those in the DataFrame
    print(rm_confounding)
    print(' '.join(sample_string))
    print(transformation)

    matching_rows = df[
        (df['Sample Names'] == ' '.join(sample_string)) &
        (df['rm_confounding'] == rm_confounding) &
        (df['Transformation'] == transformation)
        ]

    # Fill NaN values in 'rm_confounding' with 0
    df['rm_confounding'] = df['rm_confounding'].fillna('0')


    # Adjust rm_confounding parameter: treat None as 0
    rm_confounding = '0' if rm_confounding == 'None' else rm_confounding

    
    # Condition 1: Match 'Sample Names'
    condition_1 = df['Sample Names'] == ' '.join(sample_string)

    # Condition 2: Match 'rm_confounding', directly comparing after NaN replacement
    condition_2 = df['rm_confounding'] == rm_confounding

    # Condition 2: Match 'rm_confounding', directly comparing after NaN replacement
    # Ensure the DataFrame column is the same type as the comparison value
    df['rm_confounding'] = df['rm_confounding'].astype(type(rm_confounding))
    condition_2 = df['rm_confounding'] == rm_confounding
   
    # Condition 3: Match 'Transformation', case-insensitive
    condition_3 = df['Transformation'].str.lower() == transformation.lower()
    # Combine conditions to filter rows
    matching_rows = df[condition_1 & condition_2 & condition_3]
    
    # If there are any matching rows, return the File String value(s)
    if not matching_rows.empty:
        file_string_return = matching_rows['file_string'].tolist()[0]
        return file_string_return
    
    # If there are no matching rows, return None
    return None


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


class db_res:
    def __init__(self, name):
        self.name = name
        self.updn_dict = {}

    def set_updn(self, indict, updn):
        if updn == 'up':
            self.updn_dict['up'] = indict
        if updn == 'dn':
            self.updn_dict['dn'] = indict

def serve_layout():
    if flask.has_request_context():
        return url_bar_and_content_div
    return html.Div([
        url_bar_and_content_div,
        layout_index,
        layout_page1,
    ])


app.layout = serve_layout


def parse_contents(contents, filename, date):
    """
    Parses contents uploaded by the user and returns a Dash HTML Div displaying the content.
    Handles different file types (CSV, Excel).

    Parameters:
    - contents (str): The contents of the uploaded file.
    - filename (str): The name of the uploaded file.
    - date (int): The upload timestamp.

    Returns:
    - dash.development.base_component.Component: A Div containing the parsed file data.
    """
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        dash_table.DataTable(
            df.to_dict('records'),
            [{'name': i, 'id': i} for i in df.columns]
        ),

        html.Hr(),  # horizontal line

        # For debugging, display the raw contents provided by the web browser
        html.Div('Raw Content'),
        html.Pre(contents[0:200] + '...', style={
            'whiteSpace': 'pre-wrap',
            'wordBreak': 'break-all'
        })
    ])

@app.callback(
    Output('alert_import_info', 'value'),
    Output('alert_import_info_div', 'style'),
    Output('page_proceed', 'style'),
    Input('df_counts', 'data'),
    Input('df_info', 'data'),
    Input('variable_selection1', 'value'),
    Input('variable_selection2', 'value'),
    Input('variable_selection3', 'value'))
def update_alert_import(df_counts, df_info, var1, var2, var3):
    if df_counts is not None and df_info is not None and var1 is None and var2 is None and var3 is None:
        df_counts = pd.read_json(StringIO(df_counts), orient='split')
        df_info = pd.read_json(StringIO(df_info), orient='split')
        if list(df_counts.columns) != list(df_info.index):
            list1 = set(df_counts.columns)
            list2 = set(df_info.index)
            diff1 = [item for item in list1 if item not in list2]
            # Find elements in list2 not in list1
            diff2 = [item for item in list2 if item not in list1]

            return 'Columns in Count data does not match Rows in Info data \n' \
                    'Items in counts columns but not in info index: {} \n' \
                    'Items in info index but not in counts column: {}'.format(diff1, diff2), {'display': 'inline-block'}, {
                'display': 'none'}
        else:
            return 'check', {'display': 'none'}, {'display': 'none'}

    if df_counts is not None and df_info is not None and var1 is not None and var2 is not None and var3 is not None:
        return 'check', {'display': 'none'}, {'float': 'right', 'margin-right': '30px'}
    else:
        raise PreventUpdate


@app.callback(Output('checkmark_counts_div', 'style'),
              Input('df_counts', 'data'))
def update_checkmark(df_counts):
    if df_counts is None:
        raise PreventUpdate
    else:
        #return {'display': 'inline-block', 'position': 'relative', 'top': '0px', 'left': '5px'}
        return {'display': 'flex', 'alignItems': 'center', 'justifyContent': 'right'}

@app.callback(
    Output("alert-dismiss", "hide"),
    Input("alert-button", "n_clicks"),
    State("alert-dismiss", "hide"),
    #prevent_initial_call=True,
)
def alert(n_clicks, hide):
    if n_clicks is None:
        PreventUpdate
    return not hide

@app.callback(Output('checkmark_info_div', 'style'),
              Input('df_info', 'data'))
def update_checkmark(df_counts):
    if df_counts is None:
        raise PreventUpdate
    else:
        #return {'display': 'flex', 'position': 'relative', 'top': '10px', 'left': '10px'}
        return {'display': 'flex', 'alignItems': 'center', 'justifyContent': 'right'}

@app.callback(Output('df_counts', 'data'),
              Output('alert_import', 'value'),
              Output('alert_import_div', 'style'),
              Output('counts_index_store', 'data'),
              Output('counts_columns_store', 'data'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'),)
def update_output(contents, filename, date):
    if contents is None:
        raise PreventUpdate

    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)

    if filename.endswith('csv'):
        sep = ','
    elif filename.endswith('tab'):
        sep = '\t'
    elif filename.endswith('tsv'):
        sep = '\t'
    else:
        return pd.DataFrame().to_json(), 'File name extension not recongized, accepted extensions are csv, tab, tsv', \
               {'display': 'inline-block'}, [], []

    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=sep, index_col=0)
    counts_index = list(df.index)
    counts_columns = list(df.columns)

    return df.to_json(date_format='iso', orient='split'), '', {'display': 'none'}, counts_index, counts_columns


@app.callback(Output('df_info', 'data'),
              Output('info_columns_store', 'data'),
              Output('info_index_store', 'data'),
              Input('upload-data-info', 'contents'),
              State('upload-data-info', 'filename'),
              State('upload-data-info', 'last_modified'))
def update_output(contents, filename, date):
    if contents is None:
        raise PreventUpdate

    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)

    if 'csv' in filename:
        sep = ','
    elif 'tab' in filename:
        sep = '\t'
    elif 'tsv' in filename:
        sep = '\t'
    else:
        return pd.DataFrame().to_json

    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=sep, index_col=0)

    return df.to_json(date_format='iso', orient='split'), list(df.columns), list(df.index)

@app.callback(Output("content", "children"), [Input("tabs", "active_tab")])
def switch_tab(at):
    if at == "tab-1":
        return dcc.Graph(id='indicator-graphic2')
    elif at == "tab-2":
        return dcc.Graph(id='hpa_graph')
    return html.P("This shouldn't ever be displayed...")

@app.callback(Output("enrichr_content", "children"), [Input("enrichr_tabs", "active_tab")])
def switch_tab(at):
    if at == "tab-1":
        return dcc.Graph(id='Enrichr_GO_bp_up')
    elif at == "tab-2":
        return dcc.Graph(id='Enrichr_GO_bp_dn')
    elif at == "tab-3":
        return dcc.Graph(id='Enrichr_GO_cell_up')
    elif at == "tab-4":
        return dcc.Graph(id='Enrichr_GO_cell_dn')
    elif at == "tab-5":
        return dcc.Graph(id='Enrichr_GO_mf_up')
    elif at == "tab-6":
        return dcc.Graph(id='Enrichr_GO_mf_dn')
    elif at == "tab-7":
        return dcc.Graph(id='Enrichr_kegg_up')
    elif at == "tab-8":
        return dcc.Graph(id='Enrichr_kegg_dn')
    return html.P("This shouldn't ever be displayed...")


@app.callback(
    Output('dataset_name_placeholder', 'children'),
    Input('btn_save_dataset', 'n_clicks_timestamp'),
    State('selected_data', 'children'),
    State('save_dataset_name', 'value'),
    State('exclude_dropdown', 'value'))
def save_dataset(n_clicks, selected_data, dataset_name, lExclude):
    if n_clicks is None:
        raise PreventUpdate
    if dataset_name is None:
        raise PreventUpdate

    datasets = json.loads(selected_data)
    lSamples = json.loads(datasets['samples'])
    dataset_name = dataset_name
    # print([dataset_name, ' '.join(lSamples)])
    # save_dataset_to_file([dataset_name, ' '.join(lSamples)])
    # write_sample_names_to_file(lSamples, 'test')
    write_dataset(lSamples, lExclude, id_name=dataset_name)


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == "/page-1":
        return layout_page1
    elif pathname == "/page-2":
        return layout_page2
    else:
        return layout_index


def table_type(df_column):
    # Note - this only works with Pandas >= 1.0.0

    if sys.version_info < (3, 0):  # Pandas 1.0.0 does not support Python 2
        return 'any'
    if isinstance(df_column.dtype, pd.DatetimeTZDtype):
        return 'datetime',
    elif (isinstance(df_column.dtype, pd.StringDtype) or
          isinstance(df_column.dtype, pd.BooleanDtype) or
          isinstance(df_column.dtype, pd.CategoricalDtype) or
          isinstance(df_column.dtype, pd.PeriodDtype)):
        return 'text'
    elif (isinstance(df_column.dtype, pd.SparseDtype) or
          isinstance(df_column.dtype, pd.IntervalDtype) or
          isinstance(df_column.dtype, pd.Int8Dtype) or
          isinstance(df_column.dtype, pd.Int16Dtype) or
          isinstance(df_column.dtype, pd.Int32Dtype) or
          isinstance(df_column.dtype, pd.Int64Dtype)):
        return 'numeric'
    else:
        return 'any'


@app.callback([Output('exclude_list', 'children')],
              [Input('exclude_dropdown', 'value')])
def exclude_helper(exclude, verbose=0):
    if verbose == 1:
        print('!!!!!!!!exclude!!!!!!!!!!!!')
        print(exclude)
    return [json.dumps(exclude)]


@app.callback(Output('variable_selection1', 'options'),
              Input('df_info', 'data'))
def select_var1(df_info):
    if df_info is not None:
        df_info = pd.read_json(StringIO(df_info), orient='split')
        ret = [{'label': j, 'value': j} for j in df_info.columns.to_list()]
        return ret


@app.callback(Output('variable_selection1_store', 'data'),
              Input('variable_selection1', 'value'))
def select_var1_page1(variable1_value):
    if variable1_value is None:
        raise PreventUpdate
    else:
        
        return variable1_value


@app.callback(Output('variable1_selected_dropdown', 'options'),
              Input('df_info', 'data'),
              Input('variable_selection1_store', 'data'))
def select_var1(df_info, variable1_value):
    if variable1_value is None:
        raise PreventUpdate
    else:
        df_info = pd.read_json(StringIO(df_info), orient='split')
        ret = [{'label': j, 'value': j} for j in df_info[str(variable1_value)].unique()]
        filtered_data = [item for item in ret if item["label"] is not None and item["value"] is not None]
        return filtered_data


@app.callback(Output('variable_selection2', 'options'),
              Input('df_info', 'data'))
def select_var2(df_info):
    if df_info is not None:
        df_info = pd.read_json(StringIO(df_info), orient='split')
        ret = [{'label': j, 'value': j} for j in df_info.columns.to_list()]
        return ret


@app.callback(Output('variable_selection2_store', 'data'),
              Input('variable_selection2', 'value'))
def select_var2_page2(variable2_value):
    if variable2_value is None:
        raise PreventUpdate
    else:
        
        return variable2_value


@app.callback(Output('variable2_selected_dropdown', 'options'),
              Input('df_info', 'data'),
              Input('variable_selection2_store', 'data'))
def select_var2(df_info, variable2_value):
    if variable2_value is None:
        raise PreventUpdate
    else:
        df_info = pd.read_json(StringIO(df_info), orient='split')
        ret = [{'label': j, 'value': j} for j in df_info[str(variable2_value)].unique()]
        filtered_data = [item for item in ret if item["label"] is not None and item["value"] is not None]
        return filtered_data


@app.callback(Output('variable_selection3', 'options'),
              Input('df_info', 'data'))
def select_var3(df_info):
    if df_info is not None:
        df_info = pd.read_json(StringIO(df_info), orient='split')
        ret = [{'label': j, 'value': j} for j in df_info.columns.to_list()]
        return ret


@app.callback(Output('variable3_selected_dropdown', 'options'),
              Input('df_info', 'data'),
              Input('variable_selection3_store', 'data'))
def select_var1(df_info, variable3_value):
    if variable3_value is None:
        raise PreventUpdate
    else:
        df_info = pd.read_json(StringIO(df_info), orient='split')
        ret = [{'label': j, 'value': j} for j in df_info[str(variable3_value)].unique()]
        filtered_data = [item for item in ret if item["label"] is not None and item["value"] is not None]
        return filtered_data


@app.callback(Output('variable_selection3_store', 'data'),
              Input('variable_selection3', 'value'))
def select_var3_page1(variable3_value):
    if variable3_value is None:
        raise PreventUpdate
    else:
        
        return variable3_value


@app.callback(Output('varaible_selection_div', 'style'),
              Input('df_info', 'data'))
def select_var_div(df_info):
    if df_info is None:
        raise PreventUpdate
    else:
        return {'display': 'inline-block'}


@app.callback(Output('exclude_dropdown', 'options'),
              Output('rm_confounding', 'options'),
              Output('meta_dropdown', 'options'),
              Output('gene_list', 'options'),
              Output('meta_dropdown_groupby', 'options'),
              # Input('df_counts', 'data'),
              # Input('df_info', 'data'))
              Input('counts_columns_store', 'data'),
              Input('counts_index_store', 'data'),
              Input('info_columns_store', 'data'))
def exclude_dropdown(counts_columns, counts_index, info_columns):
    if counts_columns is not None and counts_index is not None and info_columns is not None:

        exclude_options = [{'label': j, 'value': j} for j in counts_columns]
        gene_list = [{'label': j, 'value': j} for j in counts_index]
        df_info_columns = [{'label': j, 'value': j} for j in info_columns]

        return exclude_options, df_info_columns, df_info_columns, gene_list, df_info_columns
    else:
        raise PreventUpdate


@app.callback(Output('exclude_dropdown', 'value'),
              Output('variable1_selected_dropdown', 'value'),
              Output('variable2_selected_dropdown', 'value'),
              Output('variable3_selected_dropdown', 'value'),
              Input('datasets', 'value'),
              # State('exclude_list', 'children'),
              State('df_info', 'data'),
              State('variable_selection1_store', 'data'),
              State('variable_selection2_store', 'data'),
              State('variable_selection3_store', 'data'),
              PreventUpdate=True)
def select_helper(dataset, df_info, variable1, variable2, variable3):
    if dataset != 'New':
        df_info = pd.read_json(StringIO(df_info), orient='split')
        lSamples, lExclude = load_dataset(dataset)
        lSamples = lSamples.split(' ')
        lExclude = lExclude.split(' ')
        df_info_temp = df_info.loc[lSamples]
        var1_uniq = df_info_temp[variable1].unique()
        var2_uniq = df_info_temp[variable2].unique()
        var3_uniq = df_info_temp[variable3].unique()

        return lExclude, var1_uniq, var2_uniq, var3_uniq

    else:
        return [], [], [], []


@app.callback(
    Output('export_plot_clicked', 'children'),
    Input('btn_export_plot', 'n_clicks_timestamp'),
    State('intermediate-DEtable', 'children'),
    State('intermediate-table', 'children'))
def export_plot(n_clicks, indata, prefixes):
    if n_clicks is None:
        raise PreventUpdate
    else:
        if indata:

            datasets = json.loads(indata)
            df_degenes = pd.read_json(StringIO(datasets['de_table']), orient='split')
            # dTranslate = dict(df_symbol.loc[df_degenes.index]['hgnc_symbol'])
            
            df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]

            # file_string = json.loads('file_string')
            file_string = datasets['file_string']
            # prefixes = json.loads(prefixes)
            # lTotal = json.loads(prefixes['prefixes'])
            # sTotal = '_'.join(lTotal)
            # sTotal = sTotal + '_for_DEplot.tab'

            outdir = os.path.join('data', 'generated')
            de_file_for_plot = file_string + '_for_DEplot.tab'
            de_file_for_plot_path = os.path.join(outdir, de_file_for_plot)

            df_degenes.to_csv(de_file_for_plot_path, sep='\t')
            #volcano_file = file_string + '_volcano.R'
            volcano_file = file_string + '_volcano.R'
            volcano_file_path = os.path.join('data', 'scripts', volcano_file)
            volcano_template_path = os.path.join('data', 'templates', 'volcano_template.R')

            # Read in the file
            with open(volcano_template_path, 'r') as file:
                filedata = file.read()

            # Replace the target string
            #filedata = filedata.replace('infile_holder', de_file_for_plot_path)
            filedata = filedata.replace('infile_holder', de_file_for_plot_path)
            
            # Write the file out again
            with open(volcano_file_path, 'w') as file:
                file.write(filedata)
            return [f'Created plot file {volcano_file}']

        else:
            print('Run DE analysis first')
            return ['Can not create plot (run the analysis)']


# Function to only update select data table
@app.callback(
    [Output('selected_data', 'children'),
     Output('alert_selection', 'hide')],
    [Input('exclude_dropdown', 'value'),
     Input('transformation', 'value'),
     Input('variable1_selected_dropdown', 'value'),
     Input('variable2_selected_dropdown', 'value'),
     Input('variable3_selected_dropdown', 'value'),
     Input('variable_selection1_store', 'data'),
     Input('variable_selection2_store', 'data'),
     Input('variable_selection3_store', 'data')],
    [State('df_info', 'data')],
    prevent_initial_call=True)
def select_info(lExclude, transformation, variable1_dropdown, variable2_dropdown, variable3_dropdown, variable1,
                variable2, variable3, df_info):
    # variable_selection1_store = column_name_variable1 = tissue
    # variable1 = LV RV etc

    if all(var is not None for var in [variable1_dropdown, variable2_dropdown, variable3_dropdown, df_info]):
        df_info = pd.read_json(StringIO(df_info), orient='split')

        transformation = 'None' if transformation == 'Sizefactor normalization' else transformation

        df_info_temp = df_info.loc[df_info[variable1].isin(variable1_dropdown), ]
        df_info_temp = df_info_temp.loc[df_info_temp[variable2].isin(variable2_dropdown), ]
        df_info_temp = df_info_temp.loc[df_info_temp[variable3].isin(variable3_dropdown), ]

        # lExclude = lExclude.split()
        # print(lExclude)

        if len(lExclude) > 0:
            # df_info_temp = df_info_temp.loc[~df_info_temp['id_tissue'].isin(lExclude), ]
            print(lExclude)
            df_info_temp = df_info_temp[~df_info_temp.index.isin(lExclude)]

        # df_info_temp.index = df_info_temp[[variable1, variable1]].apply(lambda x: '_'.join(x), axis=1)
        # df_info_temp['full_name'] = [x+'_'+df_info_temp.loc[x, variable1] for x in df_info_temp.index]

        datasets = {variable1: json.dumps(df_info_temp[variable1].unique().tolist()),
                    'transformation': transformation,
                    variable2: json.dumps(df_info_temp[variable2].unique().tolist()),
                    variable3: json.dumps(df_info_temp[variable3].unique().tolist()),
                    'exclude': json.dumps(lExclude),
                    'samples': json.dumps(df_info_temp.index.to_list()),
                    # 'full_name': json.dumps(df_info_temp['full_name'].to_list()),
                    'empty': '0'}

        return json.dumps(datasets), ''

    else:
        empty = {'empty': '1'}
        return json.dumps(empty), 'testtesttest'


def create_de_table_comp(data, output):
    return html.Div(
        [
            dt.DataTable(
                id=output,
                data=data.to_dict('records'),
                columns=[{'id': c, 'name': c} for c in data.columns],
                # style_as_list_view=True,
                selected_rows=[],
                style_cell={'padding': '5px',
                            'whiteSpace': 'no-wrap',
                            'overflow': 'hidden',
                            'textOverflow': 'ellipsis',
                            'maxWidth': 0,
                            'maxHeight': 30,
                            'height': 30,
                            'textAlign': 'left'},
                style_header={
                    'backgroundColor': 'white',
                    'fontWeight': 'bold',
                    'color': 'black'
                },
                style_table={
                    'maxHeight': '600px',
                    'overflowY': 'scroll'}

            ),
        ], className="row", style={'margin-top': '35',
                                   'margin-left': '15',
                                   'border': '1px solid #C6CCD5'}
    )


@app.callback(
    [Output('de_table_comparison', 'children'),
     Output('table-box', 'children')],
    [Input('add_de_table', 'n_clicks_timestamp'),
     Input('clear_de_table', 'n_clicks_timestamp')],
    [State('DE-table', 'data'),
     State('de_table_comparison', 'children')]
)
def add_de_set(btn_add: int, btn_clear: int, de_table: List[dict],
               detable_comp: str) -> Tuple[str, Any]:
    if btn_add is None:
        raise PreventUpdate

    if btn_clear is not None and btn_clear > btn_add:
        return json.dumps(
            {'de_table_comp': pd.DataFrame({'hgnc': []}).to_json(orient='split', date_format='iso')}), None

    de_table = pd.DataFrame.from_dict(de_table) if de_table else pd.DataFrame()

    if de_table.empty:
        print('DE table empty')
        return json.dumps({'de_table_comp': de_table.to_json(orient='split', date_format='iso')}), None

    if detable_comp is not None:
        datasets2 = json.loads(detable_comp)
        detable_comp = pd.read_json(StringIO(datasets2['de_table_comp']), orient='split')

    de_table.index = de_table['Ensembl']
    de_table = de_table.drop(['Ensembl'], axis=1)
    de_table_sorted = de_table.sort_values(by=['log2FoldChange'], ascending=True)
    de_table_sorted['rank'] = range(1, len(de_table_sorted) + 1)

    print(de_table_sorted)
    # print(detable_comp)

    intermediate_data = json.loads(intermediate_data)

    if detable_comp is not None:
        # if comparison table is already >2
        detable_comp = detable_comp.sort_values(by=['log2FoldChange'], ascending=True)
        detable_comp = range(1, len(detable_comp) + 1)
        detable_comp = detable_comp.reindex(de_table_sorted.index)

    else:
        intermediate_data = intermediate_data.sort_values(by=['log2FoldChange'], ascending=True)
        intermediate_data = range(1, len(intermediate_data) + 1)
        intermediate_data = intermediate_data.reindex(de_table_sorted.index)

    print(de_table_sorted, intermediate_data)

    result = pd.concat([de_table_sorted, intermediate_data], axis=1)
    print(result)

    dataset = {'de_table_comp': result.to_json(orient='split', date_format='iso')}
    data = create_de_table_comp(result, 'de_table_comp')

    return json.dumps(dataset), data


background_pipe = {
    'width': '15px',
    'height': '40px',
    'backgroundColor': '#4CB5F5',
    'z-index': '1',
    'left': '20%',
    'position': 'relative',
}


##main table: intermediate-table
@app.callback(
    [Output('intermediate-table', 'children'),
     Output('alert_main_table', 'hide'),
     Output('submit-done', 'children'),
     Output('select_to_explore_block', 'style'),
     Output('sample_to_pca_block', 'style'),
     Output('pca_to_de_block', 'style')],
    Input('btn-selected_data_submit', 'n_clicks'),
    State('selected_data', 'children'),
    State('rm_confounding', 'value'),
    State('full_text', 'value'),
    State('force_run', 'value'),
    State('df_counts', 'data'),
    State('df_info', 'data'),
    State('variable_selection1_store', 'data'),
    State('variable_selection2_store', 'data'),
    State('variable_selection3_store', 'data'),
    prevent_initial_call=True)
def log2_table_update(n_clicks, selected_data, rm_confounding, fulltext, force_run, df_counts, df_info,
                      variable_selection1, variable_selection2, variable_selection3):
    if n_clicks == 0:
        raise PreventUpdate
    else:

        if df_counts is not None:
            df_counts = pd.read_json(StringIO(df_counts), orient='split')
            df_info = pd.read_json(StringIO(df_info), orient='split')
            out_folder = os.path.join('data','generated')
            if not os.path.isdir(out_folder):
                cmd = 'mkdir %s' % out_folder
                os.system(cmd)

            datasets = json.loads(selected_data)
            empty = json.loads(datasets['empty'])

            if empty != '0':
                var1_dropdown = json.loads(datasets[variable_selection1])
                var2_dropdown = json.loads(datasets[variable_selection2])
                var3_dropdown = json.loads(datasets[variable_selection3])
                exclude = json.loads(datasets['exclude'])
                transformation = datasets['transformation']
                lSamples = json.loads(datasets['samples'])

                if len(lSamples) == 0:

                    df_info_temp = df_info.loc[
                        df_info.loc[df_info[variable_selection1].isin(var1_dropdown), ].index]

                    df_info_temp = df_info_temp.loc[
                        df_info_temp.loc[df_info_temp[variable_selection2].isin(var2_dropdown), ].index]

                    df_info_temp = df_info_temp.loc[df_info_temp[variable_selection3].isin(var3_dropdown), ]

                    if exclude:
                        for sample_excl in exclude:
                            if sample_excl in df_info_temp.index:
                                df_info_temp = df_info_temp.drop(sample_excl)
                else:
                    print(lSamples)
                    # loading samples from datasets.txt
                    df_info_temp = df_info.loc[lSamples,]

                if fulltext:
                    df_counts_raw = df_counts[df_info_temp.index]
                    # df_info_temp[['id_tissue', variable_selection2]].apply(lambda x: '_'.join(x), axis=1)
                    df_info_temp.index = df_info_temp['samples']
                    # Change names in count files
                    df_counts_raw.columns = df_info_temp.index
                else:
                    df_counts_raw = df_counts[df_info_temp.index]

                if rm_confounding is None:
                    rm_confounding = 'None'

                file_string = search_session(lSamples, transformation, rm_confounding)
                print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                print(file_string, 'file_string')
                print("before write",file_string)
                new_run = False

                if file_string is None:
                    new_run = True
                    file_string = generate_random_string()
                    write_session_to_file(list(df_info_temp.index), rm_confounding, transformation, file_string)
                else:
                    print(file_string, 'is not None')
                
                name_counts_for_pca = os.path.join('data', 'generated', f'{file_string}_counts.tab')
                
                name_meta_for_pca = os.path.join('data', 'generated', f'{file_string}_meta.tab')

                # Save the dataframes to the respective files
                df_counts_raw.to_csv(name_counts_for_pca, sep='\t')
                df_info_temp.to_csv(name_meta_for_pca, sep='\t')

                # Construct the output path for normalized data
                name_out = os.path.join('data', 'generated', f'{file_string}_normalized.tab')

                # Construct the path for the R script in a cross-platform way
                r_script_path = os.path.join('functions', 'normalize_vsd_rlog_removebatch2.R')

                # Format the command to run the R script, using the constructed paths
                cmd = f'Rscript {r_script_path} {name_counts_for_pca} {name_meta_for_pca} {rm_confounding} {name_out} {transformation}'

                if force_run:
                        print(cmd)
                        os.system(cmd)
                else:
                    if new_run:
                        print(cmd)
                        os.system(cmd)
                    else:
                        print('Loading file: %s' % name_out)

                df_counts_temp_norm = pd.read_csv(name_out, sep='\t', index_col=0)

                datasets = {'counts_norm': df_counts_temp_norm.to_json(orient='split', date_format='iso'),
                            'transformation': transformation,
                            'meta': df_info_temp.to_json(orient='split', date_format='iso'),
                            'counts_raw': df_counts_raw.to_json(orient='split', date_format='iso'),
                            'counts_raw_file_name': json.dumps(name_counts_for_pca),
                            'perf_file': json.dumps(name_out),
                            'file_string': json.dumps(file_string)}
                # 'prefix_sub': json.dumps(sPrefix_tissue_batch_type)}

            alert_mess = 'Data loaded successfully.'

            return json.dumps(datasets), alert_mess, '', background_pipe, background_pipe, background_pipe

        else:
            raise PreventUpdate


@app.callback(
    Output('barplot', 'figure'),
    Input('intermediate-table', 'children'),
    Input('meta_dropdown_groupby', 'value'),
    State('variable_selection1_store', 'data'),
    State('variable_selection2_store', 'data'),
    State('variable_selection3_store', 'data'))
def update_output1(indata, meta_dropdown_groupby, variable1, variable2, variable3):
    if indata is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        df_meta_temp = pd.read_json(StringIO(datasets['meta']), orient='split')
        ltraces = []

        # Variable3 == batch variable

        for x in df_meta_temp[meta_dropdown_groupby].unique():
            df_temp = df_meta_temp.loc[df_meta_temp[meta_dropdown_groupby] == x, ]
            df_temp = df_temp.groupby(variable3).count()
            trace = go.Bar(x=df_temp.index, y=df_temp['id_tissue'], name=x)
            ltraces.append(trace)

        return {
            'data': ltraces,
            'layout': go.Layout(
            )
        }
    

@app.callback(
    [Output('volcanoplot', 'figure')],
    [Input('volcanoplot-input', 'value'),
     Input('intermediate-DEtable', 'children'),
     Input('pvalue', 'children')],
    [State('volcano_xaxis', 'value'),
     State('volcano_yaxis', 'value')])
def generate_volcano(effects, indata, psig, xaxis, yaxis):
    if indata is None:
        raise PreventUpdate
    else:
        if psig is None:
            psig = 0.05

        datasets = json.loads(indata)
        df_volcano = pd.read_json(StringIO(datasets['de_table']), orient='split')
        df_volcano['-log10(p)'] = df_volcano['padj'].apply(float)
        df_volcano['pvalue'] = df_volcano['pvalue'].apply(float)

        xlabel = xaxis

        if '-log10(p)' in yaxis:
            genomewideline_value = -math.log10(float(psig))
            highlight = True
            logp = True

        else:
            genomewideline_value = False
            highlight = False
            logp = False

        radiode = datasets['DE_type']
        title_ = datasets['file_string']
        # pval = 'padj'

        pval = yaxis
        effect_size = xaxis
        print('HIGHLIGHT ', highlight)
        dashVolcano = dashbio.VolcanoPlot(
            dataframe=df_volcano,
            genomewideline_value=genomewideline_value,
            effect_size_line=effects,
            effect_size=effect_size,
            p=pval,
            gene='Ensembl',
            logp=logp,
            snp='hgnc',
            #title=title_,
            xlabel=xlabel,
            highlight=highlight)
        return [dashVolcano]


@app.callback(
    Output('biplot_text_radio', 'options'),
    [Input('biplot_radio', 'value')])
def set_biplot_options(selected):
    try:
        if selected[0] == 'biplot':
            return [{'label': k, 'value': k} for k in ['None', 'Text']]
        else:
            # If 'biplot' is not selected[0], return empty options
            return []
    except IndexError:
        # Handle the case where 'selected' is empty or not as expected
        return []
    except Exception as e:
        print(f"Unexpected error: {e}")
        return []

@app.callback(
    Output('hpa_graph', 'figure'),
    [Input('gene_list', 'value'),
     Input('radio_symbol', 'value'),
     Input('input-2_hgnc', 'value')
     ])
def update_hpa(lgenes, input_symbol, lgenes_hgnc):
    if lgenes_hgnc:
        dTranslate = dict(df_symbol_sym_ind.loc[lgenes_hgnc]['ensembl_gene_id'])
        lgenes_hgnc = [dTranslate[x] for x in lgenes_hgnc]
        lgenes = lgenes + lgenes_hgnc

    lgenes = [x for x in lgenes if x in df_hpa.index]
    hpa_table = df_hpa.loc[lgenes,]

    traces = []
    if input_symbol == 'Symbol':
        hpa_table.index = hpa_table['Gene name']
        dTranslate = dict()
        dTranslate = dict(df_symbol.loc[lgenes]['hgnc_symbol'])
        lgenes = [dTranslate.get(x, x) for x in lgenes]

    for gene in lgenes:
        traces.append(go.Bar(name=gene, x=hpa_table.loc[gene,]['Tissue'], y=hpa_table.loc[gene,]['NX']))

    return {'data': traces,
            'layout': go.Layout(title='', autosize=True, boxmode='group',
                                margin={"l": 200, "b": 100, "r": 200}, xaxis={"showticklabels": True},
                                yaxis={"title": "NX"})}


@app.callback(
    Output('indicator-graphic2', 'figure'),
    [Input('intermediate-table', 'children'),
     Input('intermediate-DEtable', 'children'),
     Input('gene_list', 'value'),
     Input('radio_symbol', 'value'),
     Input('radio-grouping', 'value'),
     Input('input-2_hgnc', 'value'),
     State('variable_selection1_store', 'data'),
     State('variable_selection2_store', 'data'),
     State('variable_selection3_store', 'data'),
     ])
def update_output1(indata, indata_de, lgenes, input_symbol, radio_grouping, lgenes_hgnc, var1, var2, var3):
    if indata is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        df_meta_temp = pd.read_json(StringIO(datasets['meta']), orient='split')
        df_counts_temp = pd.read_json(StringIO(datasets['counts_norm']), orient='split')
        sig_gene = 0
        try:
            datasets_de = json.loads(indata_de)
            df_degenes = pd.read_json(StringIO(datasets_de['de_table']), orient='split')
            df_degenes = df_degenes.loc[df_degenes['padj'] <= 0.05,]
            sig_gene = 1

        except:
            pass

        dTranslate = dict()

        if lgenes_hgnc:
            dTranslate = dict(df_symbol_sym_ind.loc[lgenes_hgnc]['ensembl_gene_id'])
            lgenes_hgnc = [dTranslate[x] for x in lgenes_hgnc]
            lgenes = lgenes + lgenes_hgnc

        if input_symbol == 'Symbol':
            dTranslate = dict(df_symbol.loc[lgenes]['hgnc_symbol'])

        traces = []
        lgrps = []
        lSamples = []

        if radio_grouping == 'var1':
            cond = var1
        elif radio_grouping == 'var2':
            cond = var2
        elif radio_grouping == 'var3':
            cond = var3
        else:
            print('error radio_grouping', radio_grouping, var1, var2, var3)
        lgrouping = df_meta_temp[cond].unique()

        '''
        if radio_grouping == 'tissue':
            lgrouping = df_meta_temp['tissue'].unique()
            cond = 'tissue'
        elif radio_grouping == 'type':
            lgrouping = df_meta_temp['type'].unique()
            cond = 'type'
        elif radio_grouping == 'SeqTag':
            lgrouping = df_meta_temp['SeqTag'].unique()
            cond = 'SeqTag'
        '''

        print(lgrouping)
        for grp in lgrouping:
            lgrps = lgrps + df_meta_temp.loc[df_meta_temp[cond] == grp,][cond].tolist()
            lSamples.append(df_meta_temp.loc[df_meta_temp[cond] == grp,].index.tolist())

        lSamples = [item for sublist in lSamples for item in sublist]

        for g in lgenes:
            try:
                signif = ''
                df_gene_t = df_counts_temp.loc[g][lSamples]
                # Mark significant genes
                if sig_gene:
                    if g in df_degenes.index:
                        signif = '*'

                traces.append(go.Box(x=lgrps, y=df_gene_t, text=df_gene_t.index,
                                     name=dTranslate.get(g, g) + signif, marker={"size": 4}, boxpoints='all',
                                     pointpos=0,
                                     jitter=0.3))
            except:
                pass
        return {
            'data': traces,
            'layout': go.Layout(title='', autosize=True, boxmode='group',
                                margin={"l": 200, "b": 100, "r": 200}, xaxis={"showticklabels": True, },
                                yaxis={"title": "log2"})}


@app.callback(
    Output('clicked', 'children'),
    [Input('pca-graph', 'clickData')])
def clicked_out(clickData):
    try:
        print(clickData['points'][0]['text'])
        return clickData['points'][0]['text']
    except:
        pass


@app.callback(
    Output('pca_and_barplot', 'figure'),
    Output("pca-graph", "figure"),
    Output('pca_comp_explained', 'figure'),
    Output('pca_correlation', 'figure'),
    Input('intermediate-table', 'children'),
    Input('meta_dropdown_groupby', 'value'),
    # Input('input-gene', 'value'),
    Input('sample_names_toggle', 'on'),
    Input('number_of_genes', 'value'),
    Input('biplot_radio', 'value'),
    Input('biplot_text_radio', 'value'),
    Input('meta_dropdown', 'value'))
def update_pca_and_barplot(indata, meta_dropdown_groupby, sample_names_toggle, number_of_genes, biplot_radio,
                           biplot_text_radio, dropdown):
    if indata is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        df_meta_temp = pd.read_json(StringIO(datasets['meta']), orient='split')
        df_counts_pca = pd.read_json(StringIO(datasets['counts_norm']), orient='split')

        if number_of_genes:
            df_counts_pca = df_counts_pca.loc[
                df_counts_pca.var(axis=1).sort_values(ascending=False).iloc[0:int(number_of_genes), ].index,]
            df_counts_pca.to_csv('~/Dropbox/dash/pca_counts.tab', sep='\t')

        print(sample_names_toggle, 'sample_names_toggle')
        pca = PCA(n_components=3)
        principalComponents = pca.fit_transform(df_counts_pca.T)
        principalDf = pd.DataFrame(data=principalComponents, columns=['PCA1', 'PCA2', 'PCA3'])
        principalDf.index = df_counts_pca.columns
        # print(principalDf)
        principalDf[dropdown] = df_meta_temp.loc[principalDf.index,][dropdown]
        traces = []
        ltraces3d = []
        lPCA_plot = []

        df_comp = pd.DataFrame(pca.components_.T)
        df_comp.index = df_counts_pca.index
        df_top = df_comp.abs().sum(axis=1).sort_values(ascending=False).iloc[0:int(4), ]
        lTop4_genes = df_top.index.tolist()

        lDropdown = df_meta_temp[dropdown].unique()

        if isinstance(lDropdown[0], (int, float)):

            if sample_names_toggle:
                mode_ = 'markers+text'
            else:
                mode_ = 'markers'  #

            traces.append(go.Scatter(x=principalDf['PCA1'], y=principalDf['PCA2'],
                                     marker={'size': 14, 'colorbar': {'title': '--'},
                                             'color': df_meta_temp[dropdown]},
                                     text=principalDf.index, mode=mode_, textposition='top center', showlegend=False))

        else:

            for grp in df_meta_temp[dropdown].unique():
                principalDf_temp = principalDf[principalDf[dropdown] == grp]

                if sample_names_toggle:
                    traces.append(
                        go.Scatter(x=principalDf_temp['PCA1'], y=principalDf_temp['PCA2'], marker={'size': 14},
                                   text=principalDf_temp.index, mode='markers+text', name=grp,
                                   textposition='top center'))
                    ltraces3d.append(go.Scatter3d(
                        x=principalDf_temp['PCA1'], y=principalDf_temp['PCA2'], z=principalDf_temp['PCA3'],
                        name=grp,
                        text=principalDf_temp.index, mode='markers+text', marker={'size': 8, 'opacity': 0.8,
                                                                                  }))

                else:
                    traces.append(
                        go.Scatter(x=principalDf_temp['PCA1'], y=principalDf_temp['PCA2'], marker={'size': 14},
                                   text=principalDf_temp.index, mode='markers', name=grp))
                    ltraces3d.append(go.Scatter3d(
                        x=principalDf_temp['PCA1'], y=principalDf_temp['PCA2'], z=principalDf_temp['PCA3'],
                        name=grp,
                        text=principalDf_temp.index, mode='markers', marker={'size': 8, 'opacity': 0.8,
                                                                             }))

        trace_var = [go.Bar(x=principalDf.columns, y=pca.explained_variance_ratio_)]

        if biplot_radio is not None:
            if len(biplot_radio) > 0:
                df_components = pd.DataFrame(pca.components_.T)
                df_components.index = df_counts_pca.index
                df_components.columns = principalDf.columns[:3]
                dTranslate = dict(df_symbol.loc[lTop4_genes]['hgnc_symbol'])
                if biplot_text_radio is not None:
                    if biplot_text_radio == 'Text':
                        mode = 'text+lines'
                    else:
                        mode = 'lines'
                else:
                    mode = 'lines'

                for gene in lTop4_genes:
                    x = [0, df_components.loc[gene, 'PCA1'] * 100]
                    y = [0, df_components.loc[gene, 'PCA2'] * 100]
                    z = [0, df_components.loc[gene, 'PCA3'] * 100]
                    hgnc_name = dTranslate.get(gene, gene)

                    if isinstance(hgnc_name, list):
                        hgnc_name = hgnc_name[1]

                    ltraces3d.append(
                        go.Scatter3d(x=x, y=y, z=z, mode=mode, marker_size=40, name=hgnc_name,
                                     text=dTranslate.get(gene, gene)))

        np.seterr(invalid='ignore')

        df_meta_for_corr = df_meta_temp.dropna(axis='columns')
        principalDf['PCA1'] = pd.to_numeric(principalDf['PCA1'], errors='coerce')
        #for column in df_meta_for_corr.columns:
        #    df_meta_for_corr[column] = pd.to_numeric(df_meta_for_corr[column], errors='coerce')

        df_meta_for_corr_numeric = df_meta_for_corr.select_dtypes(include=[np.number])

        #correlation_matrix = df_meta_for_corr.corrwith(principalDf['PCA1'], method='pearson')

        correlation_matrix = df_meta_for_corr_numeric.corrwith(principalDf['PCA1'], method='pearson')
        
        lPCA_data = []
        lPCA_data.append({
            'z': correlation_matrix,
            'y': ['PCA1' for x in correlation_matrix.index.tolist()],
            'x': correlation_matrix.index.tolist(),
            'reversescale': 'true',
            'colorscale': [[0, 'white'], [1, 'blue']],
            'type': 'heatmap',
        })

        #correlation_matrix2 = df_meta_temp.corr(method='pearson')
        #correlation_matrix2 = correlation_matrix2.fillna(0)
        #lPCA_data2 = []
        #lPCA_data2.append({
        #    'z': correlation_matrix2,
        #    'y': correlation_matrix2.index.tolist(),
        #    'x': correlation_matrix2.columns.tolist(),
        #    'reversescale': 'true',
        #    'colorscale': [[0, 'white'], [1, 'blue']],
        #    'type': 'heatmap',
        #})

        width_cor = correlation_matrix.shape[0] * 10
        
        return {'data': traces,
                'layout': go.Layout(title='Expression values', autosize=True, boxmode='group',
                                    margin={"l": 200, "b": 100, "r": 200},
                                    xaxis={'title': 'PCA1', "showticklabels": False},
                                    yaxis={"title": "PCA2"})}, \
               {"data": ltraces3d,
                "layout": go.Layout(
                    height=700, title="...",
                    # paper_bgcolor="#f3f3f3",
                    scene={"aspectmode": "cube", "xaxis": {"title": "PCA1 %.3f" % pca.explained_variance_ratio_[0], },
                           "yaxis": {"title": "PCA2 %.3f" % pca.explained_variance_ratio_[1], },
                           "zaxis": {"title": "PCA3 %.3f" % pca.explained_variance_ratio_[2], }},
                    clickmode='event+select'), }, \
               {'data': trace_var,
                'layout': go.Layout(title='', autosize=True, boxmode='group',
                                    margin={"l": 200, "b": 100, "r": 200}, xaxis={"showticklabels": True},
                                    yaxis={"title": "Variance explained"})}, \
               {'data': lPCA_data,
                'layout': {
                    'height': 100,
                    'width': 500,
                    'xaxis': {'side': 'top'},
                    'margin': {
                        'l': 200,
                        'r': 200,
                        'b': 150,
                        't': 100
                    }
                }
                }


@app.callback(Output('datasets', 'options'),
              Input('dataset_loader_start', 'children'))
def populate_dataset_load(invalue):
    dataset_file_path = os.path.join('data','datasets','datasets.csv')

    if os.path.exists(dataset_file_path):
        df_datasets = pd.read_csv(dataset_file_path)
        if not df_datasets.empty:
            lDatasets = df_datasets['ID Name'].tolist()
            return [{'label': j, 'value': j} for j in lDatasets]
        else:
            return []
    else:
        return []
    # options = [{'label': j, 'value': j} for j in lPapers], value = 'New'


@app.callback(
    [Output('DE-table', 'data'),
     Output('pvalue', 'children'),
     Output('number_of_degenes', 'children'),
     Output('de_to_enrichr_block', 'style'), ],
    [Input('sig_submit', 'n_clicks')],
    [State('volcanoplot-input', 'value'),
     State('intermediate-DEtable', 'children'),
     State('toggle_sig', 'value'),
     State('toggle_basemean', 'value')], )
def update_de_table(n_clicks, effects, indata, sig_value, basemean):
    if n_clicks is None:
        raise PreventUpdate
    else:
        if indata is None:
            raise PreventUpdate
        else:
            datasets = json.loads(indata)
            print('update_de_table')
            df_degenes = pd.read_json(StringIO(datasets['de_table']), orient='split')
            # print(df_degenes)
            radiode = datasets['DE_type']
            df_degenes = df_degenes.loc[df_degenes['padj'] <= float(sig_value),]
            if 'baseMean' in df_degenes.columns:
                df_degenes = df_degenes.loc[df_degenes['baseMean'] >= float(basemean),]

            overlap_genes = list(set(df_degenes.index).intersection(set(df_symbol.index)))
            dTranslate2 = dict(df_symbol.loc[overlap_genes]['wikigene_description'])
            df_degenes['wikigene_desc'] = [dTranslate2.get(x, x) for x in df_degenes.index]
            df_degenes1 = df_degenes.loc[df_degenes['log2FoldChange'] < effects[0],]
            df_degenes2 = df_degenes.loc[df_degenes['log2FoldChange'] > effects[1],]

            df_degenes = pd.concat([df_degenes1, df_degenes2])
            number_of_degenes = df_degenes.shape[0]
            # number_of_degenes = '#DE-genes: ' + str(number_of_degenes) + ' / ' + str(df_counts_combined.shape[0]), \
            #                    ' (%s | %s)'%(df_degenes1.shape[0], df_degenes2.shape[0]), ' Fold Change: %s %s'\
            #                    %(effects[0], effects[1])

            name = f"{datasets['file_string']}_{effects[0]}_{effects[1]}_de.tab"
            overlap_genes = list(set(df_degenes.index).intersection(set(df_symbol.index)))
            dTranslate = dict(df_symbol.loc[overlap_genes]['hgnc_symbol'])
            df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
            output_path = os.path.join('data', 'generated', name)
            df_degenes.to_csv(output_path, sep='\t')

            return df_degenes.to_dict('records'), sig_value, number_of_degenes, background_pipe


@app.callback(
    Output('export_placeholder', 'children'),
    Input('btn_export', 'n_clicks'),
    State('intermediate-DEtable', 'children'),
    State('intermediate-table', 'children'))
def export_de(n_clicks, indata, prefixes):
    if n_clicks is None:
        raise PreventUpdate
    else:
        if indata:

            datasets = json.loads(indata)
            df_degenes = pd.read_json(StringIO(datasets['de_table']), orient='split')
            df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
            file_string = datasets['file_string']
            de_file_for_plot = file_string + '_for_DEplot.tab'
            output_path = os.path.join('data', 'generated', de_file_for_plot)
            df_degenes.to_csv(output_path, sep='\t')
            print('Written to file: %s' % de_file_for_plot)

            return ['Written to file: %s' % de_file_for_plot]

        else:
            print('Run DE analysis first')
            return ['Nothing to export']


# DE ANALYSIS
@app.callback(
    Output('intermediate-DEtable', 'children'),
    Output('temp', 'children'),
    Output('submit-de-done', 'children'),
    Input('btn-DE', 'n_clicks'),
    State('intermediate-table', 'children'),
    State('program', 'value'),
    State('transformation', 'value'),
    State('force_run', 'on'),
    State('rowsum', 'value'),
    State('design', 'value'),
    State('reference', 'value'))
def run_DE_analysis(n_clicks, indata, program, transformation, force_run, rowsum, design, reference):
    if n_clicks is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        # df_meta_temp = pd.read_json(datasets['meta'], orient='split')
        outdir = os.path.join('data','generated')
        name_counts = json.loads(datasets['counts_raw_file_name'])
        file_string = json.loads(datasets['file_string'])
        logfile = os.path.join(outdir, file_string + '_DE.log')

        data = {
            'program': [program],
            'transformation': [transformation],
            'rowsum': [rowsum],
            'design': [design],
            'reference': [reference],
            'ID': [file_string]
        }
        df_new = pd.DataFrame(data)
        df_new['timestamp'] = datetime.now()
        print('1')
        if os.path.isfile(logfile):
            df_log = pd.read_csv(logfile)
            #df_log = df_log.append(df_new, ignore_index=True)
            df_log = pd.concat([df_log, df_new], ignore_index=True)
        else:
            # if logfile does not exist, use the new dataframe
            df_log = df_new

        df_log.to_csv(logfile, index=False)
        print('2')
        # lTotal = json.loads(datasets['prefixes'])
        # sTotal = '_'.join(lTotal)
        # file_string_program = lTotal + [program]
        # sTotalDE = sTotal + '_' + program
        name_meta = os.path.join('data', 'generated', f'{file_string}_meta.tab')
        name_out = os.path.join('data', 'generated', f'{file_string}_DE.tab')

        '''
        elif program == 'limma':
            if sTotal.count('@') > 1:
                cmd = 'Rscript ./functions/run_limma_de.R %s %s %s %s %s' % (
                name_counts, name_meta, 'batch', reference,
                name_out)
            else:
                cmd = 'Rscript ./functions/run_limma_de.R %s %s %s %s %s' % (
                name_counts, name_meta, 'none', reference,
                name_out)
        '''
        print('3')
        if not os.path.isfile(name_out):
            #print(cmd)
            #os.system(cmd)
            if program == 'DESeq2':
                r_script_path = os.path.join('functions', 'run_deseq2_v2.R')
                cmd = f'Rscript {r_script_path} {name_counts} {name_meta} {rowsum} {design} {reference} {name_out}'

            elif program == 'edgeR':
                r_script_path = os.path.join('functions', 'run_edgeR.R')
                cmd = f'Rscript {r_script_path} {name_counts} {name_meta} {design} {reference} {name_out}'

        else:
            if force_run:
                print(cmd)
                os.system(cmd)

        df_degenes = pd.read_csv(name_out, sep='\t')

        if program == 'edgeR':
            df_degenes.rename(columns={'logFC': 'log2FoldChange', 'FDR': 'padj', 'PValue': 'pvalue'}, inplace=True)
            df_degenes['Ensembl'] = df_degenes['genes']

        df_degenes['Ensembl'] = df_degenes.index
        df_degenes = df_degenes.sort_values(by=['log2FoldChange'])
        df_degenes['padj'] = df_degenes['padj'].apply(str)
        df_degenes['pvalue'] = df_degenes['pvalue'].apply(str)
        overlap_genes = list(set(df_degenes.index).intersection(set(df_symbol.index)))
        dTranslate = dict(df_symbol.loc[overlap_genes]['hgnc_symbol'])

        df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
        datasets = {'de_table': df_degenes.to_json(orient='split'), 'DE_type': program, 'file_string': file_string}
        print('4')
        print('DE done')
        return json.dumps(datasets), 'temp', ''


@app.callback(
    Output('Enrichr_GO_bp_up', 'figure'),
    Input('Enrichr_GO_bp_up_ph', 'figure')
)
def enrichr_tab_out(in1):
    return in1


@app.callback(
    Output('Enrichr_GO_bp_dn', 'figure'),
    Input('Enrichr_GO_bp_dn_ph', 'figure')
)
def enrichr_tab_out(in1):
    return in1


@app.callback(
    Output('Enrichr_GO_cell_up', 'figure'),
    Input('Enrichr_GO_cell_up_ph', 'figure')
)
def enrichr_tab_out(in1):
    return in1


@app.callback(
    Output('Enrichr_GO_cell_dn', 'figure'),
    Input('Enrichr_GO_cell_dn_ph', 'figure')
)
def enrichr_tab_out(in1):
    return in1


@app.callback(
    Output('Enrichr_GO_mf_up', 'figure'),
    Input('Enrichr_GO_mf_up_ph', 'figure')
)
def enrichr_tab_out(in1):
    return in1


@app.callback(
    Output('Enrichr_GO_mf_dn', 'figure'),
    Input('Enrichr_GO_mf_dn_ph', 'figure')
)
def enrichr_tab_out(in1):
    return in1


@app.callback(
    Output('Enrichr_kegg_up', 'figure'),
    Input('Enrichr_kegg_up_ph', 'figure')
)
def enrichr_tab_out(in1):
    return in1


@app.callback(
    Output('Enrichr_kegg_dn', 'figure'),
    Input('Enrichr_kegg_dn_ph', 'figure')
)
def enrichr_tab_out(in1):
    return in1


##Enrichr
@app.callback(
    [Output('Enrichr_GO_bp_up_ph', 'figure'),
     Output('Enrichr_GO_bp_dn_ph', 'figure'),
     Output('Enrichr_GO_cell_up_ph', 'figure'),
     Output('Enrichr_GO_cell_dn_ph', 'figure'),
     Output('Enrichr_GO_mf_up_ph', 'figure'),
     Output('Enrichr_GO_mf_dn_ph', 'figure'),
     Output('Enrichr_kegg_up_ph', 'figure'),
     Output('Enrichr_kegg_dn_ph', 'figure'),
     Output('submit_done_enrichr', 'children')],
    Input('btn-enrichr', 'n_clicks'),
    State('intermediate-DEtable', 'children'),
    State('DE-table', 'data'),
    State('pvalue', 'children'))
def enrichr_up(n_clicks, indata, indata_de, sig_value):
    if n_clicks is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        df_degenes = pd.DataFrame.from_dict(indata_de)

        radiode = datasets['DE_type']
        file_string = datasets['file_string']
        #
        l_databases = ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                       'KEGG_2016']
        lOutput = []
        lUpDn = ['up', 'dn']
        out_folder = os.path.join('data','generated','enrichr')

        prefix = os.path.join(out_folder, file_string)

        if not os.path.isdir(out_folder):
            cmd = 'mkdir %s' % out_folder
            os.system(cmd)

        df_degenes = df_degenes.loc[df_degenes['padj'] <= float(sig_value),]
        # else:
        # df_degenes = df_degenes.loc[df_degenes['pvalue'] <= float(sig_value), ]

        df_degenes_up = df_degenes.loc[df_degenes['log2FoldChange'] > 0,]
        df_degenes_dn = df_degenes.loc[df_degenes['log2FoldChange'] < 0,]

        genes_path_all = prefix + '_DE_genes.txt'
        wf = open(genes_path_all, 'w')
        for hgnc in df_degenes['hgnc']:
            if not hgnc == None:
                if isinstance(hgnc, list):
                    hgnc = hgnc[0]
                wf.write(hgnc + '\n')
        wf.close()

        for db in l_databases:
            for updn in lUpDn:
                if updn == 'up':
                    df_degenes = df_degenes_up
                else:
                    df_degenes = df_degenes_dn

                genes_path = prefix + '_%s_DE_genes.txt' % updn

                if len(df_degenes['hgnc']) > 5:

                    # if not os.path.isfile(genes_path):
                    wf = open(genes_path, 'w')
                    for hgnc in df_degenes['hgnc']:
                        if not hgnc == None:
                            if isinstance(hgnc, list):
                                hgnc = hgnc[0]
                            wf.write(hgnc + '\n')
                    wf.close()
                    fOut = prefix + '_' + updn + '_' + db

                    # Dont perform the same analysis again
                    if not os.path.isfile(fOut):
                        cmd = 'python3 ./enrichr-api/query_enrichr_py3.py %s %s %s %s' % (
                            genes_path, updn + '_' + db, db, fOut)
                        os.system(cmd)

                    df = pd.read_csv(fOut + '.txt', sep='\t', index_col=None)
                    df = df.sort_values(by='Adjusted P-value')
                    df_sig = len(df.loc[df['Adjusted P-value'] <= 0.05])

                    if df_sig > 0:
                        df_10 = df.head(df_sig)
                        df_10.index = df_10['Term']
                        df_10['-log(padj)'] = df_10['Adjusted P-value'].map(lambda a: -np.log(a))
                        df_10['Genes involved (%)'] = df_10['Overlap'].map(
                            lambda a: 100 * (int(a.split('/')[0]) / int(a.split('/')[1])))
                        trace = [go.Bar(x=df_10['Term'], y=df_10['Genes involved (%)'],
                                        marker={'color': df_10['Adjusted P-value'],
                                                'colorscale': 'Magma', 'showscale': True})]
                    else:
                        trace = [go.Bar()]

                    output = {'data': trace,
                              'layout': go.Layout(title=db + '_' + updn)
                              }
                else:
                    trace = [go.Bar()]
                    output = {'data': trace,
                              'layout': go.Layout(title=db + '_' + updn)
                              }

                lOutput.append(output)

    return lOutput[0], lOutput[1], lOutput[2], lOutput[3], lOutput[4], lOutput[5], lOutput[6], lOutput[7], ''

if __name__ == "__main__":
    # app.run_server(debug=True, host='130.238.239.158')
    ascii_banner = pyfiglet.figlet_format("RNA analysis")
    print(ascii_banner)
    app.run_server(debug=True, host='localhost')
