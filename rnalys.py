# -*- coding: utf-8 -*-
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
import io

from dash import dcc
from dash import html
from dash import Input, Output, State
from dash import dash_table
from collections import Counter
from dash.exceptions import PreventUpdate
from sklearn.decomposition import PCA
from demos import dash_reusable_components as drc

from layout_content import layout_index, layout_page1




'''
New in version 8:
    Include wgcna analysis v.0.1
    
Version 9
    added smooth quantile normalization and limma voom DE analysis

Version 10
    added cytoscape page with enricher and network analysis

Version 11
    added paired sample alternative in design (only works for tissue atm)
    
Version 13 
    add design formula
    removed paired sample

version 15
    remove cytoscape

version 16
    front-end changes
    cleaned up
    
'''


# external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

app = dash.Dash(
    __name__,
    #external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css', dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
    #external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'],
    external_stylesheets=[dbc.themes.LITERA],
    suppress_callback_exceptions=True
)

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),

    #Data shared between pages
    dcc.Store(id='df_counts'),#, storage_type='session'),
    dcc.Store(id='df_info'),# storage_type='session'),
    dcc.Store(id='variable_selection1_store'),
    dcc.Store(id='variable_selection2_store'),
    dcc.Store(id='variable_selection3_store'),

    html.Div(id='page-content')
])



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

lExclude = []
lTissue = []
PAGE_SIZE = 10

#Load Human protein atlas (consensus rna tissue)
df_hpa = pd.read_csv('data/rna_tissue_consensus.tsv', sep='\t', index_col='Gene')

#Setup for Enrichr
databases_enrichr = ['TRANSFAC_and_JASPAR_PWMs', 'OMIM_Disease', 'OMIM_Expanded', 'KEGG_2016', 'BioCarta_2016',
                     'WikiPathways_2016', 'Panther_2016', 'GO_Biological_process_2018', 'GO_Cellular_Component_2018',
                     'GO_Molecular_Function_2018', 'Human_Phenotype_Ontology', 'MGI_Mammalian_Phenotype_2017',
                     'Jensen_DISEASES']

#Create temporary df for volcanoplot init
columns = ['Ensembl', 'hgnc', 'baseMean', 'log2FoldChange', 'pvalue', 'padj']
data = [['ens', 'hgnc1', 10, 4, 0.005, 0.9995]]
df = pd.DataFrame(data, columns=columns)

layout_index = layout_index
layout_page1 = layout_page1


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

def save_dataset_to_file(lFilename_list):
    dataset_file = open('data/datasets.txt', 'a')
    dataset_file.write(lFilename_list[0]+' '+lFilename_list[1]+ "\n")
    dataset_file.close()

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

#TODO check compability between info and counts data

#@app.callback(Output('df_combability_check'), 'value')

@app.callback(
            Output('alert_import_info', 'value'),
            Output('alert_import_info_div', 'style'),
            Output('page_proceed', 'style'),
            Input('df_counts', 'data'),
            Input('df_info', 'data'),
            Input('variable_selection1','value'),
            Input('variable_selection2','value'),
            Input('variable_selection3','value'))
def update_alert_import(df_counts, df_info, var1, var2, var3):
    if df_counts is not None and df_info is not None and var1 is not None and var2 is not None and var3 is not None:
        return 'check', {'display': 'none'}, {'display': 'inline-block'}
    else:
        raise PreventUpdate

@app.callback(Output('checkmark_counts_div', 'style'),
              Input('df_counts', 'data'))
def update_checkmark(df_counts):
    if df_counts is None:
        raise PreventUpdate
    else:
        return {'display':'inline-block'}


@app.callback(Output('checkmark_info_div', 'style'),
              Input('df_info', 'data'))
def update_checkmark(df_counts):
    if df_counts is None:
        raise PreventUpdate
    else:
        return {'display':'flex'}#{'display':'inline-block'}


@app.callback(Output('df_counts', 'data'),
              Output('alert_import', 'value'),
              Output('alert_import_div', 'style'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'))
def update_output(contents, filename, date):
    if contents is None:
        raise PreventUpdate

    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)

    if filename.endswith('csv'):
        sep=','
    elif filename.endswith('tab'):
        sep = '\t'
    elif filename.endswith('tsv'):
        sep = '\t'
    else:
        return pd.DataFrame().to_json(), 'File name extension not recongized, accepted extensions are csv, tab, tsv', {'display': 'inline-block'}

    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=sep, index_col=0)

    return df.to_json(date_format='iso', orient='split'), '', {'display': 'none'}

@app.callback(Output('df_info', 'data'),
              Input('upload-data-info', 'contents'),
              State('upload-data-info', 'filename'),
              State('upload-data-info', 'last_modified'))
def update_output(contents, filename, date):
    if contents is None:
        raise PreventUpdate

    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)

    if 'csv' in filename:
        sep=','
    elif 'tab' in filename:
        sep = '\t'
    elif 'tsv' in filename:
        sep = '\t'
    else:
        return pd.DataFrame().to_json

    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=sep, index_col=0)
    return df.to_json(date_format='iso', orient='split')



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
    State('save_dataset_name', 'value'))
def save_dataset(n_clicks, selected_data, dataset_name):
    if n_clicks is None:
        raise PreventUpdate

    datasets = json.loads(selected_data)
    lSamples = json.loads(datasets['samples'])
    dataset_name = ''.join(dataset_name.split())
    print([dataset_name, ' '.join(lSamples)])
    save_dataset_to_file([dataset_name, ' '.join(lSamples)])

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == "/page-1":
        return layout_page1
    elif pathname == "/page-2":
        return layout_page2
    else:
        return layout_index

def load_dataset(paper):
    datasets_file = 'data/datasets.txt'

    if os.path.isfile(datasets_file):
        with open(datasets_file, 'r') as f:
            for line in f:
                if line.split()[0] == paper:
                    samples = line.split()[1:]
    return samples

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
              [Input('exclude', 'value')])
def exclude_helper(exclude, verbose=0):
    if verbose == 1:
        print('!!!!!!!!!!!!!!!!!!!!')
        print(exclude)
    return [json.dumps(exclude)]


#   html.Div(id='varaible_selection_div', children=[
#        html.Div([
#            html.P('Select Tissue column'),
#            dcc.Dropdown(id='varable_selection_div')
#        ]),

@app.callback(Output('variable_selection1', 'options'),
             Input('df_info', 'data'))
def select_var1(df_info):
    if df_info is not None:
        df_info = pd.read_json(df_info, orient='split')
        ret = [{'label': j, 'value': j} for j in df_info.columns.to_list()]
        return ret

@app.callback(Output('variable_selection1_store', 'data'),
             Input('variable_selection1', 'value'))
def select_var1_page1(variable1_value):
    if variable1_value is None:
        raise PreventUpdate
    else:
        print('variable1_value', variable1_value)
        return variable1_value


@app.callback(Output('variable1_selected_dropdown', 'options'),
             Input('df_info', 'data'),
             Input('variable_selection1_store', 'data'))
def select_var1(df_info, variable1_value):
    if variable1_value is None:
        raise PreventUpdate
    else:
        df_info = pd.read_json(df_info, orient='split')
        ret = [{'label': j, 'value': j} for j in df_info[str(variable1_value)].unique()]
        filtered_data = [item for item in ret if item["label"] is not None and item["value"] is not None]
        return filtered_data

@app.callback(Output('variable_selection2', 'options'),
             Input('df_info', 'data'))
def select_var2(df_info):

    if df_info is not None:
        df_info = pd.read_json(df_info, orient='split')
        ret = [{'label': j, 'value': j} for j in df_info.columns.to_list()]
        return ret

@app.callback(Output('variable_selection2_store', 'data'),
             Input('variable_selection2', 'value'))
def select_var2_page2(variable2_value):
    if variable2_value is None:
        raise PreventUpdate
    else:
        print('variable2_value', variable2_value)
        return variable2_value


@app.callback(Output('variable2_selected_dropdown', 'options'),
             Input('df_info', 'data'),
             Input('variable_selection2_store', 'data'))
def select_var2(df_info, variable2_value):
    if variable2_value is None:
        raise PreventUpdate
    else:
        df_info = pd.read_json(df_info, orient='split')
        ret = [{'label': j, 'value': j} for j in df_info[str(variable2_value)].unique()]
        filtered_data = [item for item in ret if item["label"] is not None and item["value"] is not None]
        return filtered_data


@app.callback(Output('variable_selection3', 'options'),
             Input('df_info', 'data'))
def select_var3(df_info):
    if df_info is not None:
        df_info = pd.read_json(df_info, orient='split')
        ret = [{'label': j, 'value': j} for j in df_info.columns.to_list()]
        return ret

@app.callback(Output('variable3_selected_dropdown', 'options'),
             Input('df_info', 'data'),
             Input('variable_selection3_store', 'data'))
def select_var1(df_info, variable3_value):
    if variable3_value is None:
        raise PreventUpdate
    else:
        df_info = pd.read_json(df_info, orient='split')
        ret = [{'label': j, 'value': j} for j in df_info[str(variable3_value)].unique()]
        filtered_data = [item for item in ret if item["label"] is not None and item["value"] is not None]
        return filtered_data

#variable_selection3

@app.callback(Output('variable_selection3_store', 'data'),
             Input('variable_selection3', 'value'))
def select_var3_page1(variable3_value):
    if variable3_value is None:
        raise PreventUpdate
    else:
        print('variable3_value', variable3_value)
        return variable3_value

@app.callback(Output('varaible_selection_div', 'style'),
             Input('df_info', 'data'))
def select_var_div(df_info):
    if df_info is None:
        raise PreventUpdate
    else:
        return {'display':'inline-block'}

@app.callback([Output('exclude', 'value'),
              Output('variable3_selected_dropdown', 'value'),
               Output('variable1_selected_dropdown', 'value'),
               Output('variable2_selected_dropdown', 'value')],
              [Input('paper', 'value')],
              [State('exclude_list', 'children'),
               State('df_info', 'data')])
def select_helper(paper, data, df_info,verbose=0):

    if df_info is not None:
        df_info = pd.read_json(df_info, orient='split')

        if paper == 'New':
            lSamples = []
            lExclude = []
            df_info_temp = df_info
            lBatch = []
            lTissue = []
            lType = []

        else:
            lSamples = load_dataset(paper)
            df_info_temp = df_info[df_info['id_tissue'].isin(lSamples)]
            lBatch = list(df_info_temp['SeqTag'].unique())
            lTissue = list(df_info_temp['tissue'].unique())
            lType = list(df_info_temp['type'].unique())

        if data:
            exclude_list = json.loads(data)

        df_info_temp = df_info[df_info['SeqTag'].isin(lBatch)]
        df_info_temp = df_info_temp[df_info_temp['tissue'].isin(lTissue)]
        df_info_temp = df_info_temp[df_info_temp['type'].isin(lType)]
        lExclude = list(df_info_temp[~df_info_temp['id_tissue'].isin(lSamples)].index)

        if verbose == 1:
            print('/////////////////////////////////////7')
            print(df_meta_temp)
            print(lTissue)
            print(lType)
            print(lBatch)
            print('/////////////////////////////////////7')

        return lExclude, lBatch, lTissue, lType
    else:
        raise PreventUpdate

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
            df_degenes = pd.read_json(datasets['de_table'], orient='split')
            #dTranslate = dict(df_symbol.loc[df_degenes.index]['hgnc_symbol'])
            #print(df_degenes.index)
            df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
            prefixes = json.loads(prefixes)
            lTotal = json.loads(prefixes['prefixes'])
            sTotal = '_'.join(lTotal)
            sTotal_deplot = sTotal + '_for_DEplot.tab'
            df_degenes.to_csv('data/generated/%s' %sTotal_deplot, sep='\t')
            sTotal = sTotal + '_volcano.R'

            # Read in the file
            with open('./data/templates/volcano_template.R', 'r') as file:
                filedata = file.read()

            # Replace the target string
            filedata = filedata.replace('infile_holder', sTotal_deplot)

            # Write the file out again
            with open('./data/scripts/%s' %sTotal, 'w') as file:
                file.write(filedata)
            return ['Created plot file %s' %sTotal]

        else:
            print('Run DE analysis first')
            return ['Can not create plot (run the analysis)']


#@app.callback(
#    Output('output-data-upload', 'children'),
#    Input('df_counts', 'data')
#)
#def output_from_store(stored_data):
#    df = pd.read_json(stored_data, orient='split')


#Function to only update select data table
@app.callback(
    [Output('selected_data', 'children'),
    Output('alert_selection', 'hide')],
    [Input('exclude', 'value'),
    Input('variable_selection3_store', 'data'),
    Input('variable_selection1_store', 'data'),
    Input('variable_selection2_store', 'data'),
    Input('transformation', 'value'),
    Input('variable3_selected_dropdown', 'value'),
    Input('variable1_selected_dropdown', 'value'),
    Input('variable2_selected_dropdown', 'value')],
    [State('df_info', 'data')],
    prevent_initial_call=True)
def select_info(lExclude, batch_dropdown,  tissue_dropdown, type_dropdown, transformation, variable3, variable1, variable2, df_info):
    #variable_selection1_store = column_name_variable1 = tissue
    #variable1 = LV RV etc
    if batch_dropdown is not None and tissue_dropdown is not None and type_dropdown is not None and transformation is not None:

        df_info = pd.read_json(df_info, orient='split')

        '''
        if 'LVRV' in tissue_dropdown:
            tissue_dropdown = tissue_dropdown + ['LV', 'RV']
        if 'EFSF' in tissue_dropdown:
            tissue_dropdown = tissue_dropdown + ['EF', 'SF']
        if 'HF' in type_dropdown:
            type_dropdown = type_dropdown + ['HFrEF', 'HFPEF']
        '''


        if transformation =='Sizefactor normalization':
            transformation = 'None'
        else:
            transformation = transformation

        df_info_temp = df_info.loc[df_info[batch_dropdown].isin(variable3),]
        #df_info_temp = df_info.loc[df_info['SeqTag'].isin(batch_dropdown),]
        #df_info_temp = df_info_temp.loc[df_info_temp['tissue'].isin(tissue_dropdown),]
        df_info_temp = df_info_temp.loc[df_info_temp[tissue_dropdown].isin(variable1),]
        #df_info_temp = df_info_temp.loc[df_info_temp['type'].isin(type_dropdown), ]
        df_info_temp = df_info_temp.loc[df_info_temp[type_dropdown].isin(variable2),]

        if lExclude is not None:
            df_info_temp = df_info_temp.loc[~df_info_temp['id_tissue'].isin(lExclude), ]


        datasets = {'tissues': json.dumps(df_info_temp['tissue'].unique().tolist()),
                    'transformation': transformation,
                    'types': json.dumps(df_info_temp['type'].unique().tolist()),
                    'batches': json.dumps(df_info_temp['SeqTag'].unique().tolist()),
                    'exclude': json.dumps(lExclude),
                    'samples': json.dumps(df_info_temp['id_tissue'].to_list()),
                    'empty': '0'}

        return json.dumps(datasets), ''


    else:
        empty = {'empty': '1'}
        return json.dumps(empty), 'testtesttest'

def change_to_hf(instr):
    if 'HFPEF' in instr:
        return 'HF'
    if 'HFrEF' in instr:
        return 'HF'
    else:
        return instr


def change_to_LVRV(instr):
    if 'LV' in instr:
        return 'LVRV'
    if 'RV' in instr:
        return 'LVRV'
    else:
        return instr



def change_to_EFSF(instr):
    if 'EF' in instr:
        return 'EFSF'
    if 'SF' in instr:
        return 'EFSF'
    else:
        return instr


def create_de_table_comp(data, output):
    return html.Div(
        [
            dt.DataTable(
                id=output,
                data=data.to_dict('records'),
                columns=[{'id': c, 'name': c} for c in data.columns],
                #style_as_list_view=True,
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
     State('de_table_comparison', 'children'),
     State('intermediate-table', 'children')],)
def add_de_set(btn_add, btn_clear, indata_de, detable_comp, indata):

    if btn_add is None:
        raise PreventUpdate
    else:
        #print(btn_add, btn_clear)
        if btn_clear is None or btn_add > btn_clear:
            de_table = pd.DataFrame.from_dict(indata_de)

            if de_table is None:
                print('DE table empty')
                df = pd.DataFrame()

            else:

                de_table.index = de_table['Ensembl']

                try:
                    datasets2 = json.loads(detable_comp)
                    de_table_comp = pd.read_json(datasets2['de_table_comp'], orient='split')
                except:
                    de_table_comp = None

                datasets3 = json.loads(indata)
                name = '_'.join(json.loads(datasets3['prefix_sub']))
                de_table = de_table.rename(columns={'log2FoldChange': name})

                if detable_comp is not None:
                    df = de_table[['hgnc', name]].copy()
                    df['hgnc_temp'] = df['hgnc']
                    df = df.drop(['hgnc'], axis=1)
                    df = pd.concat([de_table_comp, df], axis=1)
                    df['hgnc'].fillna(df['hgnc_temp'], inplace=True)
                    df = df.drop(['hgnc_temp'], axis=1)

                else:
                    df = de_table[['hgnc', name]].copy()

        else:
            if btn_clear > btn_add:
                df = pd.DataFrame({'hgnc': []})

    dataset = {'de_table_comp': df.to_json(orient='split', date_format='iso')}
    data = create_de_table_comp(df, 'de_table_comp')

    return json.dumps(dataset), data

##main table: intermediate-table
@app.callback(
    [Output('intermediate-table', 'children'),
     Output('alert_main_table', 'hide'),
     Output('submit-done', 'children')],
    Input('btn-selected_data_submit', 'n_clicks'),
    State('selected_data', 'children'),
    State('rm_confounding', 'value'),
    State('full_text', 'value'),
    State('force_run', 'value'),
    State('variable1_selected_dropdown', 'value'),
    State('df_counts', 'data'),
    State('df_info', 'data'),
    prevent_initial_call=True,)
def log2_table_update(n_clicks, selected_data, rm_confounding, fulltext, force_run, tissue_dropdn, df_counts, df_info):
    #print(n_clicks)
    if n_clicks is 0:
        raise PreventUpdate
    else:
        if df_counts is not None:
            df_counts = pd.read_json(df_counts, orient='split')
            df_info = pd.read_json(df_info, orient='split')
            out_folder = 'data/generated'
            if not os.path.isdir(out_folder):
                cmd = 'mkdir %s' % out_folder
                os.system(cmd)

            #print(selected_data)
            datasets = json.loads(selected_data)

            empty = json.loads(datasets['empty'])
            if empty != '0':

                tissue_dropdown = json.loads(datasets['tissues'])
                type_dropdown = json.loads(datasets['types'])
                batch_dropdown = json.loads(datasets['batches'])
                exclude = json.loads(datasets['exclude'])
                transformation = datasets['transformation']
                lSamples = json.loads(datasets['samples'])

                #print('tissue_dropdown:', tissue_dropdown)

                if 'HF' in type_dropdown:
                    df_info['type'] = df_info['type'].apply(change_to_hf)

                if 'LVRV' in tissue_dropdn:
                    df_info['tissue'] = df_info['tissue'].apply(change_to_LVRV)

                if 'EFSF' in tissue_dropdn:
                    df_info['tissue'] = df_info['tissue'].apply(change_to_EFSF)

                if len(lSamples) == 0:
                    # Fiter tissue
                    df_info_temp = df_info.loc[
                        df_info.loc[df_info['tissue'].isin(tissue_dropdown),].index]
                    # Filter type
                    df_info_temp = df_info_temp.loc[
                        df_info_temp.loc[df_info_temp['type'].isin(type_dropdown),].index]
                    # Filter batch
                    df_info_temp = df_info_temp.loc[df_info_temp['SeqTag'].isin(batch_dropdown),]

                    if exclude:
                        for sample_excl in exclude:
                            if sample_excl in df_info_temp.index:
                                df_info_temp = df_info_temp.drop(sample_excl)
                else:
                    print(lSamples)
                    #loading samples from datasets.txt
                    df_info_temp = df_info.loc[lSamples,]

                if fulltext:
                    df_counts_raw = df_counts[df_info_temp.index]
                    df_info_temp.index = df_info_temp[['id_tissue', 'type']].apply(lambda x: '_'.join(x), axis=1)
                    #Change names in count files
                    df_counts_raw.columns = df_info_temp.index
                else:
                    df_counts_raw = df_counts[df_info_temp.index]

                if exclude:
                    exclude = [x.split('SLL')[1] for x in exclude]
                else:
                    exclude = ['']
                #print(df_meta_combined['tissue'])
                lTotal = exclude + tissue_dropdown + batch_dropdown + type_dropdown
                sPrefix_tissue_batch_type = tissue_dropdown + batch_dropdown + type_dropdown+[transformation]

                name_counts_for_pca = './data/generated/' + ''.join(lTotal) + '_counts.tab'
                name_meta_for_pca = './data/generated/' + ''.join(lTotal) + '_meta.tab'

                #remember
                #dTranslate = dict(df_symbol['hgnc_symbol'])
                #df_counts_raw.index = [dTranslate.get(x, x) for x in df_counts_raw.index]


                df_counts_raw.to_csv(name_counts_for_pca, sep='\t')
                df_info_temp.to_csv(name_meta_for_pca, sep='\t')

                # Setup the design for linear model
                if rm_confounding != None:
                    performance = 'batch_'+str(rm_confounding)

                else:
                    performance = 'batch_nobatch'

                #Run R script for normalization and transformation
                print('PERFOMRAMCE', performance)
                name_out = './data/generated/' + ''.join(lTotal) + '_' + performance + '_' + transformation + '_normalized.tab'
                cmd = 'Rscript ./functions/normalize_vsd_rlog_removebatch2.R %s %s %s %s %s' % (
                    name_counts_for_pca, name_meta_for_pca, performance, name_out, transformation)

                if not os.path.isfile(name_out):
                    print(cmd)
                    os.system(cmd)
                else:
                    if force_run:
                        print(cmd)
                        os.system(cmd)
                    else:
                        print('Loading file: %s' %name_out)

                df_counts_temp_norm = pd.read_csv(name_out, sep='\t', index_col=0)

                datasets = {'counts_norm': df_counts_temp_norm.to_json(orient='split', date_format='iso'),
                            'transformation': transformation,
                            'meta': df_info_temp.to_json(orient='split', date_format='iso'),
                            'counts_raw': df_counts_raw.to_json(orient='split', date_format='iso'),
                            'counts_raw_file_name': json.dumps(name_counts_for_pca),
                            'perf_file': json.dumps(name_out),
                            'prefixes': json.dumps(lTotal),
                            'prefix_sub': json.dumps(sPrefix_tissue_batch_type)}

                outf = '_'.join(lTotal)
            alert_mess = 'test out'

            return json.dumps(datasets), alert_mess, ''

        else:
            raise PreventUpdate

@app.callback(
    Output('barplot', 'figure'),
    [Input('intermediate-table', 'children'),
     Input('radio_coloring', 'value'),])
def update_output1(indata, radio_coloring):
    if indata is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        df_meta_temp = pd.read_json(datasets['meta'], orient='split')
        ltraces = []

        if radio_coloring == 'tissue':
            for tissue in df_meta_temp['tissue'].unique():
                df_tissue = df_meta_temp.loc[df_meta_temp['tissue'] == tissue,]
                df_batch_tissue = df_tissue.groupby('SeqTag').count()
                trace = go.Bar(x=df_batch_tissue.index, y=df_batch_tissue['id_tissue'], name=tissue)
                ltraces.append(trace)


        elif radio_coloring == 'type':
            for type in df_meta_temp['type'].unique():
                df_type = df_meta_temp.loc[df_meta_temp['type'] == type,]
                df_batch_type = df_type.groupby('type').count()
                trace = go.Bar(x=df_batch_type.index, y=df_batch_type['id_tissue'], name=type)
                ltraces.append(trace)

        elif radio_coloring == 'SeqTag':
            for batch in df_meta_temp['SeqTag'].unique():
                df_batch = df_meta_temp.loc[df_meta_temp['SeqTag'] == batch,]
                df_batch = df_batch.groupby('SeqTag').count()
                trace = go.Bar(x=df_batch.index, y=df_batch['id_tissue'], name='test')
                ltraces.append(trace)

        else:
            #Prevent update on gene-level
            raise PreventUpdate

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
        df_volcano = pd.read_json(datasets['de_table'], orient='split')
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
        title_ = datasets['title']
        #pval = 'padj'

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
            title=title_,
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
    except:
        #Biplot not selected
        return []

@app.callback(
    Output('hpa_graph', 'figure'),
    [Input('input-2', 'value'),
     Input('radio_symbol', 'value'),
     Input('input-2_hgnc', 'value')
    ])
def update_hpa(lgenes, input_symbol, lgenes_hgnc):

    if lgenes_hgnc:

        dTranslate = dict(df_symbol_sym_ind.loc[lgenes_hgnc]['ensembl_gene_id'])
        lgenes_hgnc = [dTranslate[x] for x in lgenes_hgnc]
        lgenes = lgenes + lgenes_hgnc

    lgenes = [x for x in lgenes if x in df_hpa.index]
    hpa_table = df_hpa.loc[lgenes, ]

    traces = []
    if input_symbol == 'Symbol':
        hpa_table.index = hpa_table['Gene name']
        dTranslate = dict()
        dTranslate = dict(df_symbol.loc[lgenes]['hgnc_symbol'])
        lgenes = [dTranslate.get(x, x) for x in lgenes]

    for gene in lgenes:
        traces.append(go.Bar(name=gene, x=hpa_table.loc[gene, ]['Tissue'], y=hpa_table.loc[gene, ]['NX']))

    return {'data': traces,
     'layout': go.Layout(title='', autosize=True, boxmode='group',
                         margin={"l": 200, "b": 100, "r": 200}, xaxis={"showticklabels": True},
                         yaxis={"title": "NX"})}

@app.callback(
    Output('indicator-graphic2', 'figure'),
    [Input('intermediate-table', 'children'),
     Input('intermediate-DEtable', 'children'),
     Input('input-2', 'value'),
     Input('radio_symbol', 'value'),
     Input('radio-grouping', 'value'),
     Input('input-2_hgnc', 'value')
    ])
def update_output1(indata, indata_de, lgenes, input_symbol, radio_grouping, lgenes_hgnc):
    if indata is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)

        df_meta_temp = pd.read_json(datasets['meta'], orient='split')
        df_counts_temp = pd.read_json(datasets['counts_norm'], orient='split')
        sig_gene = 0
        try:
            datasets_de = json.loads(indata_de)
            df_degenes = pd.read_json(datasets_de['de_table'], orient='split')
            df_degenes = df_degenes.loc[df_degenes['padj'] <= 0.05, ]
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

        if radio_grouping == 'tissue':
            lgrouping = df_meta_temp['tissue'].unique()
            cond = 'tissue'
        elif radio_grouping == 'type':
            lgrouping = df_meta_temp['type'].unique()
            cond = 'type'
        elif radio_grouping == 'SeqTag':
            lgrouping = df_meta_temp['SeqTag'].unique()
            cond = 'SeqTag'

        print(lgrouping)
        for grp in lgrouping:
            lgrps = lgrps + df_meta_temp.loc[df_meta_temp[cond] == grp,][cond].tolist()
            lSamples.append(df_meta_temp.loc[df_meta_temp[cond] == grp,].index.tolist())

        lSamples = [item for sublist in lSamples for item in sublist]

        for g in lgenes:
            try:
                signif = ''
                df_gene_t = df_counts_temp.loc[g][lSamples]
                #Mark significant genes
                if sig_gene:
                    if g in df_degenes.index:
                        signif = '*'

                traces.append(go.Box(x=lgrps, y=df_gene_t, text=df_gene_t.index,
                                     name=dTranslate.get(g, g)+signif, marker={"size": 4}, boxpoints='all', pointpos=0,
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
    Input('radio_coloring', 'value'),
    Input('input-gene', 'value'),
    Input('radio-text', 'value'),
    Input('num_genes', 'value'),
    Input('biplot_radio', 'value'),
    Input('biplot_text_radio', 'value'),
    Input('meta_dropdown', 'value'))
def update_pca_and_barplot(indata, radio_coloring, input2, radio_text, number_of_genes, biplot_radio,
                           biplot_text_radio, dropdown):

    if indata is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        df_meta_temp = pd.read_json(datasets['meta'], orient='split')
        df_counts_pca = pd.read_json(datasets['counts_norm'], orient='split')

        if number_of_genes:
            df_counts_pca = df_counts_pca.loc[
                df_counts_pca.var(axis=1).sort_values(ascending=False).iloc[0:int(number_of_genes), ].index,]
            df_counts_pca.to_csv('~/Dropbox/dash/pca_counts.tab', sep='\t')


        pca = PCA(n_components=3)
        principalComponents = pca.fit_transform(df_counts_pca.T)
        principalDf = pd.DataFrame(data=principalComponents, columns=['PCA1', 'PCA2', 'PCA3'])
        principalDf.index = df_counts_pca.columns
        #print(principalDf)
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

            if radio_text == 'Text':
                mode_ = 'markers+text'
            else:
                mode_ = 'markers'#

            traces.append(go.Scatter(x=principalDf['PCA1'], y=principalDf['PCA2'],
                                     marker={'size': 14, 'colorbar': {'title':'--'},
                                             'color': df_meta_temp[dropdown]},
                                     text=principalDf.index, mode=mode_, textposition='top center', showlegend=False))

        else:

            for grp in df_meta_temp[dropdown].unique():
                principalDf_temp = principalDf[principalDf[dropdown] == grp]


                if radio_text == 'Text':
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
        correlation_matrix = df_meta_for_corr.corrwith(principalDf['PCA1'], method='pearson')

        lPCA_data = []
        lPCA_data.append({
            'z': correlation_matrix,
            'y': ['PCA1' for x in correlation_matrix.index.tolist()],
            'x': correlation_matrix.index.tolist(),
            'reversescale': 'true',
            'colorscale': [[0, 'white'], [1, 'blue']],
            'type': 'heatmap',
        })

        correlation_matrix2 = df_meta_temp.corr(method='pearson')
        correlation_matrix2 = correlation_matrix2.fillna(0)
        lPCA_data2 = []
        lPCA_data2.append({
            'z': correlation_matrix2,
            'y': correlation_matrix2.index.tolist(),
            'x': correlation_matrix2.columns.tolist(),
            'reversescale': 'true',
            'colorscale': [[0, 'white'], [1, 'blue']],
            'type': 'heatmap',
        })

        width_cor = correlation_matrix.shape[0]*10
        return {'data': traces,
                'layout': go.Layout(title='Expression values', autosize=True, boxmode='group',
                                    margin={"l": 200, "b": 100, "r": 200},
                                    xaxis={'title': 'PCA1', "showticklabels": False},
                                    yaxis={"title": "PCA2"})}, \
               {"data": ltraces3d,
                "layout": go.Layout(
                    height=700, title="...",
                    #paper_bgcolor="#f3f3f3",
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
                        'height': 500,
                        'width': 1500,
                        'xaxis': {'side': 'top'},
                        'margin': {
                            'l': 200,
                            'r': 200,
                            'b': 150,
                            't': 100
                        }
                    }
                }

@app.callback(
    [Output('DE-table', 'data'),
     Output('pvalue', 'children'),
     Output('number_of_degenes', 'children')],
    [Input('sig_submit', 'n_clicks')],
    [State('volcanoplot-input', 'value'),
     State('intermediate-DEtable', 'children'),
     State('toggle_sig', 'value'),
     State('toggle_basemean', 'value')],)
def update_de_table(n_clicks, effects, indata, sig_value, basemean):
    if n_clicks is None:
        raise PreventUpdate
    else:
        if indata is None:
            raise PreventUpdate
        else:
            datasets = json.loads(indata)
            print('update_de_table')
            df_degenes = pd.read_json(datasets['de_table'], orient='split')
            #print(df_degenes)
            radiode = datasets['DE_type']
            df_degenes = df_degenes.loc[df_degenes['padj'] <= float(sig_value), ]
            df_degenes = df_degenes.loc[df_degenes['baseMean'] >= float(basemean), ]
            overlap_genes = list(set(df_degenes.index).intersection(set(df_symbol.index)))
            dTranslate2 = dict(df_symbol.loc[overlap_genes]['wikigene_description'])
            df_degenes['wikigene_desc'] = [dTranslate2.get(x, x) for x in df_degenes.index]
            df_degenes1 = df_degenes.loc[df_degenes['log2FoldChange'] < effects[0], ]
            df_degenes2 = df_degenes.loc[df_degenes['log2FoldChange'] > effects[1], ]

            df_degenes = pd.concat([df_degenes1, df_degenes2])
            number_of_degenes = df_degenes.shape[0]
            number_of_degenes = '#DE-genes: ' + str(number_of_degenes) + ' / ' + str(df_counts_combined.shape[0]), \
                                ' (%s | %s)'%(df_degenes1.shape[0], df_degenes2.shape[0]), ' Fold Change: %s %s'\
                                %(effects[0], effects[1])

            name = datasets['title']+'_'+str(effects[0])+'_'+str(effects[1]) + '_de.tab'
            overlap_genes = list(set(df_degenes.index).intersection(set(df_symbol.index)))
            dTranslate = dict(df_symbol.loc[overlap_genes]['hgnc_symbol'])
            df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
            df_degenes.to_csv('data/generated/%s' %name, sep='\t')

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
            df_degenes = pd.read_json(datasets['de_table'], orient='split')
            df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
            prefixes = json.loads(prefixes)
            lTotal = json.loads(prefixes['prefixes'])
            sTotal = '_'.join(lTotal)
            sTotal = sTotal + '_for_DEplot.tab'
            df_degenes.to_csv('data/generated/%s' %sTotal, sep='\t')
            print('Written to file: %s' %sTotal)

            return ['Written to file: %s' %sTotal]

        else:
            print('Run DE analysis first')
            return ['Nothing to export']


#DE ANALYSIS
@app.callback(
    Output('intermediate-DEtable', 'children'),
    Output('temp', 'children'),
    Output('submit-de-done', 'children'),
    Input('btn-DE', 'n_clicks'),
    State('intermediate-table', 'children'),
    State('program', 'value'),
    State('transformation', 'value'),
    State('force_run', 'value'),
    State('rowsum','value'),
    State('design', 'value'),
    State('reference', 'value'))
def run_DE_analysis(n_clicks, indata, radiode, transformation, force_run, rowsum, design, reference):
    if n_clicks is None:
        raise PreventUpdate
    else:
        datasets = json.loads(indata)
        #df_meta_temp = pd.read_json(datasets['meta'], orient='split')

        name_counts = json.loads(datasets['counts_raw_file_name'])
        lTotal = json.loads(datasets['prefixes'])
        sTotal = '_'.join(lTotal)
        lTotalDE = lTotal + [radiode]
        sTotalDE = sTotal + '_' + radiode

        name_meta = './data/generated/' + ''.join(lTotal) + '_meta.tab'
        name_out = './data/generated/' + ''.join(lTotalDE) + '_' + transformation + '_DE.tab'
        print(sTotal)
        if radiode == 'DESeq2':
            # counts, meta file, transformation, rowsum, name_out
            print(design)
            cmd = 'Rscript ./functions/run_deseq2_v2.R %s %s %s %s %s %s %s' % (
                name_counts, name_meta, transformation, rowsum, design, name_out, reference)

        elif radiode == 'edgeR':
            if sTotal.count('@') > 1:
                cmd = 'Rscript ./functions/run_edgeR.R %s %s %s %s %s' % (name_counts, name_meta, 'batch', radiode_type,
                                                                          name_out)
                print('edgeR batch')
            else:
                cmd = 'Rscript ./functions/run_edgeR.R %s %s %s %s %s' % (name_counts, name_meta, 'none', radiode_type,
                                                                          name_out)
                print('edgeR no batch')

        elif radiode == 'limma':
            if sTotal.count('@') > 1:
                cmd = 'Rscript ./functions/run_limma_de.R %s %s %s %s %s' % (
                name_counts, name_meta, 'batch', radiode_type,
                name_out)
            else:
                cmd = 'Rscript ./functions/run_limma_de.R %s %s %s %s %s' % (
                name_counts, name_meta, 'none', radiode_type,
                name_out)
        else:
            print('ERROR radio_de')

        if not os.path.isfile(name_out):
            print(cmd)
            os.system(cmd)
        else:
            if force_run:
                print(cmd)
                os.system(cmd)

        df_degenes = pd.read_csv(name_out, sep='\t')

        if radiode == 'edgeR':
            df_degenes.rename(columns={'logFC': 'log2FoldChange', 'FDR': 'padj', 'PValue': 'pvalue'}, inplace=True)
            df_degenes['Ensembl'] = df_degenes['genes']

        df_degenes['Ensembl'] = df_degenes.index
        df_degenes = df_degenes.sort_values(by=['log2FoldChange'])
        df_degenes['padj'] = df_degenes['padj'].apply(str)
        df_degenes['pvalue'] = df_degenes['pvalue'].apply(str)
        overlap_genes = list(set(df_degenes.index).intersection(set(df_symbol.index)))
        dTranslate = dict(df_symbol.loc[overlap_genes]['hgnc_symbol'])

        #logging
        #if

        df_degenes['hgnc'] = [dTranslate.get(x, x) for x in df_degenes.index]
        datasets = {'de_table': df_degenes.to_json(orient='split'), 'DE_type': radiode, 'title': sTotalDE}

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
    Output('Enrichr_kegg_dn_ph', 'figure')],
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
        #df_degenes = pd.read_json(datasets['de_table'], orient='split')
        radiode = datasets['DE_type']
        sTotalDE = datasets['title']
        #print(df_degenes.head())
        l_databases = ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                      'KEGG_2016']
        lOutput = []
        lUpDn = ['up', 'dn']
        out_folder = 'data/generated/enrichr'

        prefix = os.path.join(out_folder, sTotalDE)

        if not os.path.isdir(out_folder):
            cmd = 'mkdir %s' % out_folder
            os.system(cmd)

        df_degenes = df_degenes.loc[df_degenes['padj'] <= float(sig_value), ]
        #else:
        #df_degenes = df_degenes.loc[df_degenes['pvalue'] <= float(sig_value), ]

        df_degenes_up = df_degenes.loc[df_degenes['log2FoldChange'] > 0, ]
        df_degenes_dn = df_degenes.loc[df_degenes['log2FoldChange'] < 0, ]

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

                    #if not os.path.isfile(genes_path):
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

                    print(df['Adjusted P-value'])
                    print('### df_sig')

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

    return lOutput[0], lOutput[1], lOutput[2], lOutput[3], lOutput[4], lOutput[5], lOutput[6], lOutput[7]

if __name__ == "__main__":
    #app.run_server(debug=True, host='130.238.239.158')
    ascii_banner = pyfiglet.figlet_format("RNA analysis")
    print(ascii_banner)
    app.run_server(debug=True, host='localhost')


