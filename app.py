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
import sys
import subprocess
import logging
import plotly.express as px

logging.basicConfig(level=logging.DEBUG)

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
from functions.deseq2 import calculate_size_factors, estimate_dispersion, fit_glm_nb
from functions.edgeR import calc_norm_factors, estimate_common_dispersion, estimate_tagwise_dispersion, glm_lrt


# Configuration for external stylesheets
FONT_AWESOME = "https://use.fontawesome.com/releases/v5.10.2/css/all.css"


app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=0.8, user-scalable=no"}],
    external_stylesheets=[dbc.themes.LITERA, FONT_AWESOME],
    suppress_callback_exceptions=True
)
"""
Initializes the Dash app with specific external stylesheets.
- dbc.themes.LITERA: Bootstrap theme for styling.
- FONT_AWESOME: Font Awesome for icons.
suppress_callback_exceptions is set to True for handling exceptions in callbacks.
"""

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

    matching_rows = df[
        (df['Sample Names'] == ' '.join(sample_string)) &
        (df['rm_confounding'] == rm_confounding) &
        (df['Transformation'] == transformation)
        ]

    df['rm_confounding'] = df['rm_confounding'].fillna('0')
    rm_confounding = '0' if rm_confounding == 'None' else rm_confounding
    condition_1 = df['Sample Names'] == ' '.join(sample_string)
    condition_2 = df['rm_confounding'] == rm_confounding
    df['rm_confounding'] = df['rm_confounding'].astype(type(rm_confounding))
    condition_2 = df['rm_confounding'] == rm_confounding
    condition_3 = df['Transformation'].str.lower() == transformation.lower()

    matching_rows = df[condition_1 & condition_2 & condition_3]
    if not matching_rows.empty:
        file_string_return = matching_rows['file_string'].tolist()[0]
        return file_string_return
    
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


def generate_example_data(type):
    # generate type: counts or meta data
    # Generating the RNASeq data using the negative binomial distribution

    number_samples = 40
    number_genes = 25000
    genes = ['Gene_' + str(i) for i in range(number_genes)]

    # Creating fixed sample IDs, with first half LV and second half RV
    samples = ['SLL_' + str(i) for i in range(number_samples)]
    tissues = ['LV'] * (number_samples // 2) + ['RV'] * (number_samples // 2)

    # Parameters for the negative binomial distribution differ by tissue type
    mean_expression_LV = 100  # Mean count for LV
    mean_expression_RV = 120  # Mean count for RV, different from LV to simulate variability
    dispersion_LV = 0.4  # Dispersion for LV
    dispersion_RV = 0.6  # Dispersion for RV, can be different to simulate biological variability

    # Initialize the dataframe to store the expression data
    expression_data = np.zeros((number_genes, number_samples))

    # Generate data for each tissue type
    for i, tissue in enumerate(tissues):
        if tissue == 'LV':
            mean_expression = mean_expression_LV
            dispersion = dispersion_LV
        else:
            mean_expression = mean_expression_RV
            dispersion = dispersion_RV

        size = mean_expression / dispersion
        prob = size / (size + mean_expression)
        expression_data[:, i] = np.random.negative_binomial(size, prob, size=number_genes)

    df_rnaseq = pd.DataFrame(data=expression_data, index=genes, columns=samples)
    # Select 10% of the genes to have even more exaggerated differences
    num_differential_genes = int(0.1 * number_genes)
    differential_gene_indices = np.random.choice(number_genes, num_differential_genes, replace=False)

    # Adjust parameters for these genes
    for index in differential_gene_indices:
        expression_data[index, :number_samples//2] = np.random.negative_binomial(mean_expression_LV * 2, prob, size=number_samples//2)  # Double the LV mean for these genes
        expression_data[index, number_samples//2:] = np.random.negative_binomial(mean_expression_RV * 2, prob, size=number_samples//2)  # Double the RV mean for these genes


    # Metadata DataFrame creation
    metadata = {
        'Tissue': tissues,
        'Disease': np.random.choice(['Disease1', 'Disease2'], size=number_samples),
        'Batch': np.random.choice(['Batch1', 'Batch2', 'Batch3'], size=number_samples),
        'BMI': np.random.uniform(18.5, 30.0, size=number_samples),
        'Diabetes': np.random.choice(['Yes', 'No'], size=number_samples)
    }
    df_metadata = pd.DataFrame(metadata, index=samples)
        
    if type == 'counts':
        df = df_rnaseq
    else:
        df = df_metadata

    return df

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

    #Success, show proceed button
    if df_counts is not None and df_info is not None and var1 is not None and var2 is not None and var3 is not None:
        return 'check', {'display': 'none'}, {'display': 'inline-block', 'margin-right': '1px'}
    else:
        raise PreventUpdate

@app.callback(Output('checkmark_counts_div', 'style'),
              Input('df_counts', 'data'))
def update_checkmark(df_counts):
    if df_counts is None:
        raise PreventUpdate
    else:
        return {'display': 'flex', 'alignItems': 'center', 'justifyContent': 'right'}

@app.callback(
    Output("alert-dismiss", "hide"),
    Input("alert-button", "n_clicks"),
    State("alert-dismiss", "hide"),
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
        return {'display': 'flex', 'alignItems': 'center', 'justifyContent': 'right'}

@app.callback(Output('df_counts', 'data'),
              Output('alert_import', 'value'),
              Output('alert_import_div', 'style'),
              Output('counts_index_store', 'data'),
              Output('counts_columns_store', 'data'),
              Input('upload-data', 'contents'),
              Input('generate-example-data', 'n_clicks'),
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'),)
def update_counts_data(contents, example_data_btn, filename, date):

    if example_data_btn is not None:
        df = generate_example_data('counts')
        counts_index = list(df.index)
        counts_columns = list(df.columns)
        
    else:
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
              Output('alert_import_info_div2', 'style'),
              Input('upload-data-info', 'contents'),
              Input('generate-example-data', 'n_clicks'),
              State('upload-data-info', 'filename'),
              State('upload-data-info', 'last_modified'))
def update_info_data(contents, example_data_btn, filename, date):

    if example_data_btn is not None:
        df = generate_example_data('meta_data')
    else:
        if contents is None:
            raise PreventUpdate
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)

        if 'csv' in filename:
            sep = ','
        elif 'tab' in filename or 'tsv' in filename:
            sep = '\t'
        else:
            return pd.DataFrame().to_json(), [], [], {'display': 'inline-block'}
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=sep, index_col=0)

    return df.to_json(date_format='iso', orient='split'), list(df.columns), list(df.index), {'display': 'none'}

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

    if sys.version_info < (3, 0):
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

#Populate the variable_selection on the index page
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

        if pd.isna(lExclude):
            lExclude = []
        else:
            lExclude = lExclude.split(' ') if isinstance(lExclude, str) else []
        
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
            df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
            file_string = datasets['file_string']

            outdir = os.path.join('data', 'generated')
            de_file_for_plot = file_string + '_for_DEplot.tab'
            de_file_for_plot_path = os.path.join(outdir, de_file_for_plot)

            df_degenes.to_csv(de_file_for_plot_path, sep='\t')
            volcano_file = file_string + '_volcano.R'
            volcano_file_path = os.path.join('data', 'scripts', volcano_file)
            volcano_template_path = os.path.join('data', 'templates', 'volcano_template.R')

            # Read in the file
            with open(volcano_template_path, 'r') as file:
                filedata = file.read()

            # Replace the target string
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

        if len(lExclude) > 0:
            df_info_temp = df_info_temp[~df_info_temp.index.isin(lExclude)]

        datasets = {variable1: json.dumps(df_info_temp[variable1].unique().tolist()),
                    'transformation': transformation,
                    variable2: json.dumps(df_info_temp[variable2].unique().tolist()),
                    variable3: json.dumps(df_info_temp[variable3].unique().tolist()),
                    'exclude': json.dumps(lExclude),
                    'samples': json.dumps(df_info_temp.index.to_list()),
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

    result = pd.concat([de_table_sorted, intermediate_data], axis=1)
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
    [
        Output('intermediate-table', 'children'),
        Output('alert_main_table', 'is_open'),
        #Output('submit-done', 'children'),
        #Output('select_to_explore_block', 'style'),
        #Output('sample_to_pca_block', 'style'),
        #Output('pca_to_de_block', 'style'),
    ],
    [Input('btn-selected_data_submit', 'n_clicks')],
    [State('selected_data', 'children'),
     State('rm_confounding', 'value'),
     State('full_text', 'value'),
     State('force_run', 'value'),
     State('df_counts', 'data'),
     State('df_info', 'data'),
     State('variable_selection1_store', 'data'),
     State('variable_selection2_store', 'data'),
     State('variable_selection3_store', 'data')],
    prevent_initial_call=True
)
def table_update(n_clicks, selected_data, rm_confounding, fulltext, force_run, df_counts, df_info,
                      variable_selection1, variable_selection2, variable_selection3):
    if n_clicks == 0:
        raise PreventUpdate

    # Perform the long-running task
    if df_counts is not None:
        df_counts = pd.read_json(StringIO(df_counts), orient='split')
        df_info = pd.read_json(StringIO(df_info), orient='split')
        out_folder = os.path.join('data', 'generated')
        if not os.path.isdir(out_folder):
            cmd = f'mkdir {out_folder}'
            try:
                subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            except subprocess.CalledProcessError as e:
                print("Error executing command:", e)

        datasets = json.loads(selected_data)
        empty = json.loads(datasets['empty'])
        if empty != '0':
            var1_dropdown = json.loads(datasets[variable_selection1])
            var2_dropdown = json.loads(datasets[variable_selection2])
            var3_dropdown = json.loads(datasets[variable_selection3])
            exclude = json.loads(datasets['exclude'])
            transformation = datasets['transformation']
            lSamples = json.loads(datasets['samples'])

            if transformation is None:
                return json.dumps({}), True

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
                df_info_temp = df_info.loc[lSamples,]

            if fulltext:
                df_counts_raw = df_counts[df_info_temp.index]
                df_info_temp.index = df_info_temp['samples']
                df_counts_raw.columns = df_info_temp.index
            else:
                df_counts_raw = df_counts[df_info_temp.index]

            if rm_confounding is None:
                rm_confounding = 'None'

            file_string = search_session(lSamples, transformation, rm_confounding)

            new_run = False
            if file_string is None:
                new_run = True
                file_string = generate_random_string()
                write_session_to_file(list(df_info_temp.index), rm_confounding, transformation, file_string)
           
            name_counts_for_pca = os.path.join('data', 'generated', f'{file_string}_counts.tab')
            name_meta_for_pca = os.path.join('data', 'generated', f'{file_string}_meta.tab')

            df_counts_raw.to_csv(name_counts_for_pca, sep='\t')
            df_info_temp.to_csv(name_meta_for_pca, sep='\t')

            name_out = os.path.join('data', 'generated', f'{file_string}_normalized.tab')
            r_script_path = os.path.join('functions', 'data_transformation.R')

            cmd = f'Rscript {r_script_path} {name_counts_for_pca} {name_meta_for_pca} {rm_confounding} {name_out} {transformation}'

            print('forcerun', force_run)

            if force_run or new_run:
                print(cmd)
                try:
                    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                except subprocess.CalledProcessError as e:
                    print("Error executing command:", e)
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

        alert_mess = 'Data loaded successfully.'

        #block_style = {'display': 'block'}

        return json.dumps(datasets), False
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

        if meta_dropdown_groupby is None:
            meta_dropdown_groupby = variable3
        
        #barplot by batch
        for x in df_meta_temp[meta_dropdown_groupby].unique():
            df_temp = df_meta_temp.loc[df_meta_temp[meta_dropdown_groupby] == x, ]
            df_temp = df_temp.groupby(variable3).count()
            trace = go.Bar(x=df_temp.index, y=df_temp.iloc[:,1], name=x)
            ltraces.append(trace)

        return {
            'data': ltraces,
            'layout': go.Layout(height=600,
                                    margin={'l':200,
                                            'r':200,
                                            'b':150,
                                            't':100
                                            }
                                    )
            
        }
    

@app.callback(
    Output('ma_plot', 'figure'),
    Input('intermediate-DEtable', 'children'),
)
def generate_maplot(indata):
    if indata is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        df_matable = pd.read_json(StringIO(datasets['ma_table']), orient='split')

        ma_plot_fig = px.scatter(
            df_matable, 
            x='mean_norm_counts', 
            y='log2FoldChange',
            labels={
                "mean_norm_counts": "Mean of Normalized Counts (A)",
                "log2FoldChange": "Log2 Fold Change (M)"
            },
            title="MA-Plot"
        )
        return ma_plot_fig

    
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

        pval = yaxis
        effect_size = xaxis
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

    #only works for imported data containing Ensembl or hgnc id
    try:
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

    except:
        traces = []

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
            #if DE anaysis has been conducted
            datasets_de = json.loads(indata_de)
            df_degenes = pd.read_json(StringIO(datasets_de['de_table']), orient='split')
            df_degenes = df_degenes.loc[df_degenes['padj'] <= 0.05,]
            sig_gene = 1
        except:
            pass
        
        #Only works for Ensembl or HGNC id
        try:
            #If the counts data have Ensembl gene names
            dTranslate = dict()
            if lgenes_hgnc:
                dTranslate = dict(df_symbol_sym_ind.loc[lgenes_hgnc]['ensembl_gene_id'])
                lgenes_hgnc = [dTranslate[x] for x in lgenes_hgnc]
                lgenes = lgenes + lgenes_hgnc

            if input_symbol == 'Symbol':
                dTranslate = dict(df_symbol.loc[lgenes]['hgnc_symbol'])
        except:
            pass

        traces = []
        lgrps = []
        lSamples = []

        print('RADIO:::::::: ',radio_grouping)
        #radio_grouping = radio_grouping[0]
                
        grouping_map = {
            'var1': var1,
            'var2': var2,
            'var3': var3
        }

        # Retrieve the appropriate condition based on radio_grouping
        cond = grouping_map.get(radio_grouping)

        # Error handling if radio_grouping does not match any key
        if cond is None:
            print('error radio_grouping', radio_grouping, var1, var2, var3)

        lgrouping = df_meta_temp[cond].unique()

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
        return clickData['points'][0]['text']
    except:
        pass


@app.callback(
    Output('pca_and_barplot', 'figure'),
    Output("pca-graph", "figure"),
    Output('pca_comp_explained', 'figure'),
    Output('pca_correlation', 'figure'),
    Input('intermediate-table', 'children'),
    Input('sample_names_toggle', 'on'),
    Input('number_of_genes', 'value'),
    Input('biplot_radio', 'value'),
    Input('biplot_text_radio', 'value'),
    Input('meta_dropdown', 'value'),
    State('variable_selection1_store', 'data'),
    State('variable_selection3_store', 'data')
    )
def update_pca_and_barplot(indata, sample_names_toggle, number_of_genes, biplot_radio,
                           biplot_text_radio, dropdown, variable1, variable3):
    if indata is None:
        raise PreventUpdate
    else:

        if dropdown is None:
            dropdown = variable1

        datasets = json.loads(indata)
        df_meta_temp = pd.read_json(StringIO(datasets['meta']), orient='split')
        df_counts_pca = pd.read_json(StringIO(datasets['counts_norm']), orient='split')

        if number_of_genes:
            df_counts_pca = df_counts_pca.loc[
                df_counts_pca.var(axis=1).sort_values(ascending=False).iloc[0:int(number_of_genes), ].index,]
            df_counts_pca.to_csv('~/Dropbox/dash/pca_counts.tab', sep='\t')

        pca = PCA(n_components=3)
        principalComponents = pca.fit_transform(df_counts_pca.T)
        principalDf = pd.DataFrame(data=principalComponents, columns=['PCA1', 'PCA2', 'PCA3'])
        principalDf.index = df_counts_pca.columns
        principalDf[dropdown] = df_meta_temp.loc[principalDf.index,][dropdown]
        traces = []
        ltraces3d = []
        lPCA_plot = []

        df_comp = pd.DataFrame(pca.components_.T)
        df_comp.index = df_counts_pca.index
        df_top = df_comp.abs().sum(axis=1).sort_values(ascending=False).iloc[0:int(4), ]
        top_four_genes = df_top.index.tolist()

        lDropdown = df_meta_temp[dropdown].unique()

        #if isinstance(lDropdown[0], (int, float)):

        #    if sample_names_toggle:
        #        mode_ = 'markers+text'
        #    else:
        #        mode_ = 'markers'  #

        #    traces.append(go.Scatter(x=principalDf['PCA1'], y=principalDf['PCA2'],
        #                             marker={'size': 14, 'colorbar': {'title': '--'},
        #                                     'color': df_meta_temp[dropdown]},
        #                             text=principalDf.index, mode=mode_, textposition='top center', showlegend=False))

        #else:

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
                try:
                    dTranslate = dict(df_symbol.loc[top_four_genes]['hgnc_symbol'])
                except:
                    pass
                if biplot_text_radio is not None:
                    if biplot_text_radio == 'Text':
                        mode = 'text+lines'
                    else:
                        mode = 'lines'
                else:
                    mode = 'lines'

                for gene in top_four_genes:
                    x = [0, df_components.loc[gene, 'PCA1'] * 100]
                    y = [0, df_components.loc[gene, 'PCA2'] * 100]
                    z = [0, df_components.loc[gene, 'PCA3'] * 100]
                    try:
                        hgnc_name = dTranslate.get(gene, gene)

                        if isinstance(hgnc_name, list):
                            gene = hgnc_name[1]
                    except:
                        pass
                        

                    ltraces3d.append(
                        go.Scatter3d(x=x, y=y, z=z, mode=mode, marker_size=40, name=gene,
                                     text=dTranslate.get(gene, gene)))

        np.seterr(invalid='ignore')

        df_meta_for_corr = df_meta_temp.dropna(axis='columns')
        principalDf['PCA1'] = pd.to_numeric(principalDf['PCA1'], errors='coerce')
        df_meta_for_corr_numeric = df_meta_for_corr.select_dtypes(include=[np.number])
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

        width_cor = correlation_matrix.shape[0] * 10
        
        return {'data': traces,
                'layout': go.Layout(height=600, boxmode='group',
                                    margin={"l": 200, "b": 100, "r": 200},
                                    xaxis={'title': 'PCA1', "showticklabels": False},
                                    yaxis={"title": "PCA2"})}, \
               {"data": ltraces3d,
                "layout": go.Layout(
                    height=600,

                    scene={"aspectmode": "cube", "xaxis": {"title": "PCA1 %.3f" % pca.explained_variance_ratio_[0], },
                           "yaxis": {"title": "PCA2 %.3f" % pca.explained_variance_ratio_[1], },
                           "zaxis": {"title": "PCA3 %.3f" % pca.explained_variance_ratio_[2], }},
                    clickmode='event+select'), }, \
               {'data': trace_var,
                'layout': go.Layout(title='', height=600, boxmode='group',
                                    margin={"l": 200, "b": 100, "r": 200}, xaxis={"showticklabels": True},
                                    yaxis={"title": "Variance explained"})}, \
               {'data': lPCA_data,
                'layout': go.Layout(height=600,
                                    margin={'l':200,
                                            'r':200,
                                            'b':150,
                                            't':100
                                            }
                                    )
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


@app.callback(
    [Output('DE-table', 'data'),
     Output('pvalue', 'children'),
     Output('number_of_degenes', 'children')],
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
            df_degenes = pd.read_json(StringIO(datasets['de_table']), orient='split')
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
            name = f"{datasets['file_string']}_{effects[0]}_{effects[1]}_de.tab"
            overlap_genes = list(set(df_degenes.index).intersection(set(df_symbol.index)))

            try:
                dTranslate = dict(df_symbol.loc[overlap_genes]['hgnc_symbol'])
                df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
            except:
                pass
            output_path = os.path.join('data', 'generated', name)
            df_degenes.to_csv(output_path, sep='\t')

            return df_degenes.to_dict('records'), sig_value, number_of_degenes


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
            try:
                df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
            except:
                pass
            
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
    State('force_run', 'value'),
    State('rowsum', 'value'),
    State('design', 'value'),
    State('reference', 'value'))
def run_DE_analysis(n_clicks, indata, program, transformation, force_run, rowsum, design, reference):
    if n_clicks is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
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
        
        if os.path.isfile(logfile):
            df_log = pd.read_csv(logfile)
            df_log = pd.concat([df_log, df_new], ignore_index=True)
        else:
            # if logfile does not exist, use the new dataframe
            df_log = df_new

        df_log.to_csv(logfile, index=False)

        name_meta = os.path.join('data', 'generated', f'{file_string}_meta.tab')
        name_out = os.path.join('data', 'generated', f'{file_string}_DE.tab')

        name_ma_table = None

        if program == 'DESeq2':
            name_ma_table =  os.path.join('data', 'generated', f'{file_string}_maplot.tab')
            r_script_path = os.path.join('functions', 'run_deseq2.R')
            cmd = f'Rscript {r_script_path} {name_counts} {name_meta} {rowsum} {design} {reference} {name_out} {name_ma_table}'
            

        elif program == 'edgeR':
            r_script_path = os.path.join('functions', 'run_edgeR.R')
            cmd = f'Rscript {r_script_path} {name_counts} {name_meta} {design} {reference} {name_out}'
        
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

        if os.path.isfile(name_out) and not force_run:
            print("Output file already exists. Use 'force_run' to override.")
        else:
            try:
                # Run the subprocess command
                result = subprocess.run(cmd, check=True, text=True, capture_output=True)
                print("DE analysis completed successfully.")
                print("Output:", result.stdout)
            except subprocess.CalledProcessError as e:
                print("Failed to run edgeR analysis.")
                print("Error:", e.stderr)

        df_degenes = pd.read_csv(name_out, sep='\t')

        if program == 'edgeR':
            df_degenes.rename(columns={'logFC': 'log2FoldChange', 'FDR': 'padj', 'PValue': 'pvalue'}, inplace=True)
            df_degenes['Ensembl'] = df_degenes['genes']

        df_degenes['Ensembl'] = df_degenes.index
        df_degenes = df_degenes.sort_values(by=['log2FoldChange'])
        df_degenes['padj'] = df_degenes['padj'].apply(str)
        df_degenes['pvalue'] = df_degenes['pvalue'].apply(str)
        overlap_genes = list(set(df_degenes.index).intersection(set(df_symbol.index)))
        try:
            dTranslate = dict(df_symbol.loc[overlap_genes]['hgnc_symbol'])
            df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
        except:
            pass

        if name_ma_table:
            ma_table = pd.read_csv(name_ma_table, sep='\t')
        else:
            ma_table = pd.DataFrame()

        datasets = {'de_table': df_degenes.to_json(orient='split'), 'DE_type': program, 'file_string': file_string, 'ma_table': ma_table.to_json(orient='split')}
        
        print('DE done')
        return json.dumps(datasets), 'temp', ''


#For future removal of Rscript with interal deseq2 algorithms
'''
# DE ANALYSIS
@app.callback(
    Output('intermediate-DEtable', 'children'),
    Output('temp', 'children'),
    Output('submit-de-done', 'children'),
    Input('btn-DE', 'n_clicks'),
    State('intermediate-table', 'children'),
    State('program', 'value'),
    State('transformation', 'value'),
    State('force_run', 'value'),
    State('rowsum', 'value'),
    State('design', 'value'),
    State('reference', 'value'))
def run_DE_analysis(n_clicks, indata, program, transformation, force_run, rowsum, design, reference):
    if n_clicks is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        outdir = os.path.join('data', 'generated')
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

        if os.path.isfile(logfile):
            df_log = pd.read_csv(logfile)
            df_log = pd.concat([df_log, df_new], ignore_index=True)
        else:
            df_log = df_new

        df_log.to_csv(logfile, index=False)

        # Load counts and metadata
        name_meta = os.path.join('data', 'generated', f'{file_string}_meta.tab')
        name_out = os.path.join('data', 'generated', f'{file_string}_DE.tab')

        counts = pd.read_csv(name_counts, sep='\t', index_col=0)
        metadata = pd.read_csv(name_meta, sep='\t', index_col=0)

        # Remove rows with low counts if rowsum parameter is set
        if rowsum is not None and rowsum != '':
            rowsum_threshold = int(rowsum)
            counts = counts[counts.sum(axis=1) >= rowsum_threshold]

        # Prepare the design matrix
        design_factors = design.replace('~', '').split('+')
        design_factors = [factor.strip() for factor in design_factors]
        design_matrix = metadata[design_factors]
        design_matrix.index = metadata.index  # Ensure the index aligns with counts columns

        # Convert categorical variables to appropriate format
        for col in design_matrix.columns:
            if design_matrix[col].dtype == object or str(design_matrix[col].dtype).startswith('category'):
                design_matrix[col] = pd.Categorical(design_matrix[col])

        # Ensure counts columns match design matrix index
        counts = counts.loc[:, design_matrix.index]

        # Prepare the contrast or coefficient to test
        coef_name = design_factors[-1]  # Testing the last factor in the design

        # Handle reference level if specified
        if reference and ':' in reference:
            factor, ref_level = reference.split(':')
            if factor in design_matrix.columns:
                design_matrix[factor] = pd.Categorical(design_matrix[factor])
                design_matrix[factor].cat.set_categories([ref_level] + [x for x in design_matrix[factor].cat.categories if x != ref_level], inplace=True)
            else:
                print(f"Reference factor {factor} not in design matrix columns")

        # Run the DE analysis using the appropriate function
        if program == 'DESeq2':
            # Normalize counts
            size_factors = calculate_size_factors(counts)
            normalized_counts = counts.div(size_factors, axis=1)

            # Estimate dispersion
            dispersions = estimate_dispersion(counts, design_matrix, size_factors)

            # Fit the model and get results
            results_df = fit_glm_nb(counts, design_matrix, size_factors, dispersions, coef_name)

            # Adjust p-values
            results_df['padj'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
            results_df.rename(columns={'p_value': 'pvalue'}, inplace=True)
            df_degenes = results_df

        elif program == 'edgeR':
            # Normalize counts
            norm_factors = calc_norm_factors(counts)

            # Estimate dispersion
            groups = design_matrix[coef_name].values
            common_dispersion = estimate_common_dispersion(counts, groups, norm_factors)
            tagwise_dispersion = estimate_tagwise_dispersion(counts, groups, norm_factors, common_dispersion)

            # Fit the model and get results
            results_df = glm_lrt(counts, design_matrix, norm_factors, tagwise_dispersion, coef_name)

            # Adjust p-values
            results_df['padj'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
            results_df.rename(columns={'logFC': 'log2FoldChange', 'p_value': 'pvalue'}, inplace=True)
            df_degenes = results_df

        else:
            print("Unknown program specified.")
            return None, None, None

        # Save the results to name_out
        df_degenes.to_csv(name_out, sep='\t')

        # The rest of the code remains the same
        df_degenes['Ensembl'] = df_degenes.index
        df_degenes = df_degenes.sort_values(by=['log2FoldChange'])
        df_degenes['padj'] = df_degenes['padj'].apply(str)
        df_degenes['pvalue'] = df_degenes['pvalue'].apply(str)
        overlap_genes = list(set(df_degenes.index).intersection(set(df_symbol.index)))
        try:
            dTranslate = dict(df_symbol.loc[overlap_genes]['hgnc_symbol'])
            df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
        except:
            pass

        # Since we didn't compute MA plot data, we'll set ma_table to an empty DataFrame
        ma_table = pd.DataFrame()

        datasets = {'de_table': df_degenes.to_json(orient='split'), 'DE_type': program, 'file_string': file_string, 'ma_table': ma_table.to_json(orient='split')}

        print('DE analysis completed.')
        return json.dumps(datasets), 'temp', ''


'''

conditions = [
    ("GO_bp_up", "GO_bp_up_ph"),
    ("GO_bp_dn", "GO_bp_dn_ph"),
    ("GO_cell_up", "GO_cell_up_ph"),
    ("GO_cell_dn", "GO_cell_dn_ph"),
    ("GO_mf_up", "GO_mf_up_ph"),
    ("GO_mf_dn", "GO_mf_dn_ph"),
    ("kegg_up", "kegg_up_ph"),
    ("kegg_dn", "kegg_dn_ph")
]

# Function that simply returns the input as output
def redirect_input(in1):
    return in1

def create_callback(output_id, input_id):
    @app.callback(
        Output(output_id, 'figure'),
        Input(input_id, 'figure')
    )
    def redirect_input(in1):
        return in1

# Generating callbacks dynamically
for condition, placeholder in conditions:
    create_callback(f'Enrichr_{condition}', f'Enrichr_{condition}_ph')


##Enrichr
@app.callback(
    [Output(f'Enrichr_{db}_{state}_ph', 'figure') 
     for db in ("GO_bp", "GO_cell", "GO_mf", "kegg") 
     for state in ("up", "dn")],
    Input('btn-enrichr', 'n_clicks'),
    State('intermediate-DEtable', 'children'),
    State('DE-table', 'data'),
    State('pvalue', 'children')
)
def enrichr_up(n_clicks, indata, indata_de, sig_value):
    if n_clicks is None:
        raise PreventUpdate

    datasets = json.loads(indata)
    de_genes_df = pd.DataFrame.from_dict(indata_de)

    up_down_files = split_genes_by_expression(de_genes_df, sig_value)
    results = run_enrichr_analysis(datasets, up_down_files)

    return results

def split_genes_by_expression(de_genes_df, sig_value):
    de_genes_df = de_genes_df[de_genes_df['padj'] <= float(sig_value)]
    gene_splits = {
        'up': de_genes_df[de_genes_df['log2FoldChange'] > 0],
        'dn': de_genes_df[de_genes_df['log2FoldChange'] < 0]
    }
    return gene_splits

def run_enrichr_analysis(datasets, gene_splits):
    radiode = datasets['DE_type']
    file_string = datasets['file_string']
    databases = ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018', 'KEGG_2016']
    out_folder = os.path.join('data', 'generated', 'enrichr')
    os.makedirs(out_folder, exist_ok=True)

    output_figures = []
    for db in databases:
        for state, genes_df in gene_splits.items():
            if len(genes_df) > 5:
                file_path = write_genes_to_file(genes_df, file_string, state, out_folder)
                result_figure = perform_enrichr_query(file_path, db, state, file_string, out_folder)
                output_figures.append(result_figure)
            else:
                output_figures.append(default_figure(db, state))

    return output_figures

def write_genes_to_file(genes_df, file_string, state, out_folder):
    file_path = os.path.join(out_folder, f'{file_string}_{state}_DE_genes.txt')
    with open(file_path, 'w') as wf:
        for hgnc in genes_df['hgnc']:
            if isinstance(hgnc, list):
                hgnc = hgnc[0]
            wf.write(hgnc + '\n')
    return file_path

def perform_enrichr_query(file_path, db, state, file_string, out_folder):
    output_path = os.path.join(out_folder, f'{file_string}_{state}_{db}')
    if not os.path.isfile(output_path):
        enrichr_path = os.path.join('enrichr-api', 'query_enrichr_py3.py')
        cmd = f'{sys.executable} {enrichr_path} {file_path} {state}_{db} {db} {output_path}'
        run_subprocess(cmd)
    return create_figure(output_path, db, state)

def run_subprocess(cmd):
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print("Error executing command:", e)

def create_figure(output_path, db, state):
    df = pd.read_csv(output_path + '.txt', sep='\t')
    df['Term'] = df['Term'].str.replace(r' \(GO:\d+\)$', '', regex=True)
    df = df.sort_values(by='Adjusted P-value')
    significant_terms = df[df['Adjusted P-value'] <= 0.05]
    if len(significant_terms) > 0:
        return plot_significant_terms(significant_terms, db, state)
    return default_figure(db, state)

def plot_significant_terms(df, db, state):
    df['-log(padj)'] = -np.log(df['Adjusted P-value'])
    df['Genes involved (%)'] = df['Overlap'].apply(lambda x: 100 * (int(x.split('/')[0]) / int(x.split('/')[1])))
    trace = go.Bar(x=df['Term'], y=df['Genes involved (%)'], marker={'color': df['-log(padj)'], 'colorscale': 'Magma', 'showscale': True})
    layout = go.Layout(title=f'{db}_{state}', xaxis={'automargin': True}, height=650)
    return {'data': [trace], 'layout': layout}

def default_figure(db, state):
    layout = go.Layout(title=f'{db}_{state} - No significant results', xaxis={'automargin': True}, height=650)
    return {'data': [go.Bar()], 'layout': layout}

@app.callback(
    Output('export_enrichr_plot', 'children'),
    Input('btn_export_enrichr_plot', 'n_clicks_timestamp'),
    State('intermediate-DEtable', 'children'),
    State('intermediate-table', 'children'))
def export_enrichr_plot(n_clicks, indata, prefixes):
    if n_clicks is None:
        raise PreventUpdate
    else:
        if indata:
            datasets = json.loads(indata)        
            file_string = datasets['file_string']
            outdir = os.path.join('data', 'scripts')
            enrichr_file = file_string + '_enrichr.R'
            enrichr_file_path = os.path.join('data', 'scripts', enrichr_file)
            enrichr_template_path = os.path.join('data', 'templates', 'plot_enrichr.py')

            # Read in the file
            with open(enrichr_template_path, 'r') as file:
                filedata = file.read()

            # Replace the target string
            filedata = filedata.replace('_fileid_', file_string)
            
            # Write the file out again
            with open(enrichr_file_path, 'w') as file:
                file.write(filedata)
            return [f'Created enrichr plot file {enrichr_file}']

        else:
            print('Run DE Enrichr first')
            return ['Can not create plot (run the analysis)']

if __name__ == "__main__":
    ascii_banner = pyfiglet.figlet_format("\\\ Rnalys")
    print(ascii_banner)
    #app.run_server(debug=True, host='localhost')
    app.run_server(debug=True, dev_tools_ui=True, dev_tools_props_check=True)
