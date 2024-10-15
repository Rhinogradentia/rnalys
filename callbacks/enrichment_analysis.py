from dash.dependencies import Input, Output
from app import app
import plotly.graph_objs as go

# Callback for switching between different Enrichr tabs
@app.callback(
    Output('enrichr_content', 'children'),
    Input('enrichr_tabs', 'active_tab')
)
def switch_enrichr_tab(active_tab):
    if active_tab == "tab-1":
        return go.Figure(data=[])
    elif active_tab == "tab-2":
        return go.Figure(data=[])
    return "This tab does not exist."

