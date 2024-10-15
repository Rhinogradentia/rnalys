import os
import json
import pandas as pd
import base64
import io
from dash.dependencies import Input, Output, State
from app import app
from io import StringIO

# Callback for updating the counts data from file upload or example data generation
@app.callback(
    Output('df_counts', 'data'),
    Output('alert_import', 'value'),
    Output('alert_import_div', 'style'),
    Output('counts_index_store', 'data'),
    Output('counts_columns_store', 'data'),
    Input('upload-data', 'contents'),
    Input('generate-example-data', 'n_clicks'),
    State('upload-data', 'filename'),
    State('upload-data', 'last_modified'),
)
def update_counts_data(contents, example_data_btn, filename, date):
    if example_data_btn is not None:
        df = generate_example_data('counts')
        counts_index = list(df.index)
        counts_columns = list(df.columns)
    else:
        if contents is None:
            return pd.DataFrame().to_json(), 'File not recognized', {'display': 'inline-block'}, [], []

        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        if filename.endswith('csv'):
            sep = ','
        elif filename.endswith('tsv') or filename.endswith('tab'):
            sep = '\t'
        else:
            return pd.DataFrame().to_json(), 'Invalid file extension', {'display': 'inline-block'}, [], []

        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=sep, index_col=0)
        counts_index = list(df.index)
        counts_columns = list(df.columns)

    return df.to_json(date_format='iso', orient='split'), '', {'display': 'none'}, counts_index, counts_columns


def generate_example_data(type):
    # Function to generate example data
    number_samples = 40
    number_genes = 25000
    genes = ['Gene_' + str(i) for i in range(number_genes)]
    samples = ['Sample_' + str(i) for i in range(number_samples)]

    expression_data = pd.DataFrame(data=np.random.poisson(100, size=(number_genes, number_samples)),
                                   index=genes, columns=samples)

    if type == 'counts':
        return expression_data
    return pd.DataFrame()
