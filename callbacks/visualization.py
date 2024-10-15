import plotly.graph_objs as go
import pandas as pd
import numpy as np
from dash.dependencies import Input, Output
from app import app
from io import StringIO
import json

# Example callback for updating a barplot
@app.callback(
    Output('barplot', 'figure'),
    Input('intermediate-table', 'children'),
    Input('meta_dropdown_groupby', 'value')
)
def update_barplot(indata, groupby_var):
    if indata is None:
        raise PreventUpdate
    
    datasets = json.loads(indata)
    df_meta_temp = pd.read_json(StringIO(datasets['meta']), orient='split')
    ltraces = []

    for x in df_meta_temp[groupby_var].unique():
        df_temp = df_meta_temp[df_meta_temp[groupby_var] == x]
        trace = go.Bar(x=df_temp.index, y=df_temp['value'], name=x)
        ltraces.append(trace)

    return {
        'data': ltraces,
        'layout': go.Layout(title='Bar Plot', xaxis={'title': groupby_var}, yaxis={'title': 'Value'})
    }
