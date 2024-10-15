from dash import html, dcc

# Layout for the index page
layout_index = html.Div([
    # Location component to track URL changes
    dcc.Location(id='url', refresh=False),

    # Navigation bar
    html.Nav([
        html.A("Home", href="/", className="nav-link"),
        html.A("Page 1", href="/page-1", className="nav-link"),
    ], className="nav"),

    # Main content area where page content will be rendered
    html.Div(id='page-content'),

    # Hidden stores for intermediate data
    dcc.Store(id='df_counts'),
    dcc.Store(id='df_info'),
    dcc.Store(id='variable_selection1_store'),
    dcc.Store(id='variable_selection2_store'),
    dcc.Store(id='variable_selection3_store'),

    # Placeholder for split dataframes to speed up referencing
    dcc.Store(id='counts_index_store'),
    dcc.Store(id='counts_columns_store'),
    dcc.Store(id='info_columns_store'),
    dcc.Store(id='info_index_store'),
])
