import dash
from dash import dcc, html, Input, Output
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
import pickle

# Load processed data
df_umap = pd.read_pickle('processed_data/df_umap.pkl')
df_tsne = pd.read_pickle('processed_data/df_tsne.pkl')

with open('processed_data/bubble_df_dict.pkl', 'rb') as f:
    bubble_df_dict = pickle.load(f)

with open('processed_data/cell_composition_dict.pkl', 'rb') as f:
    cell_composition_dict = pickle.load(f)

with open('processed_data/color_discrete_maps.pkl', 'rb') as f:
    color_discrete_maps = pickle.load(f)

with open('processed_data/n_clusters_options.pkl', 'rb') as f:
    n_clusters_options = pickle.load(f)

# Create initial pie chart
def create_pie_chart(cell_counts, title):
    labels = cell_counts.index
    values = cell_counts.values
    
    # Get color mapping from the UMAP plot for consistency
    unique_cell_types = df_umap['Cell Class'].unique()
    color_map = {cell_type: px.colors.qualitative.D3[i % len(px.colors.qualitative.D3)] 
                 for i, cell_type in enumerate(unique_cell_types)}
    
    # Map colors to the current labels
    colors = [color_map[label] for label in labels]
    
    # Calculate percentages
    total = sum(values)
    percentages = [v/total*100 for v in values]
    
    # Create text template that only shows labels for percentages > 1%
    texttemplate = ['%{label}<br>%{percent:.1%}' if p > 1 else '' for p in percentages]
    
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        hole=.3,
        texttemplate=texttemplate,
        textposition='outside',
        textinfo='text',
        marker_colors=colors,
        showlegend=False
    )])
    
    fig.update_layout(
        title=dict(
            text=f'<b>{title}</b>',
            x=0.5,
            xanchor='center',
            font=dict(size=20, family='Inter')
        ),
        height=400,
        margin=dict(l=20, r=20, t=50, b=20),
        font=dict(family='Inter'),
        template='simple_white'
    )
    
    return fig

# Get initial pie chart data
initial_cluster = '0'  # Start with first cluster
initial_cell_counts = cell_composition_dict[('UMAP', n_clusters_options[0])][initial_cluster]
fig_pie = create_pie_chart(initial_cell_counts, f'Cell Type Composition<br>UMAP Cluster {initial_cluster}')

# Initialize Dash app with custom index string to include Inter font
app = dash.Dash(__name__, 
    index_string='''
    <!DOCTYPE html>
    <html>
        <head>
            <link rel="preconnect" href="https://rsms.me/">
            <link rel="stylesheet" href="https://rsms.me/inter/inter.css">
            {%metas%}
            <title>{%title%}</title>
            {%favicon%}
            {%css%}
            <style>
                * {
                    font-family: 'Inter', sans-serif;
                }
                body {
                    background-color: white !important;
                    margin: 0;
                    padding: 20px;
                }
                #_dash-app-content {
                    background-color: white !important;
                }
            </style>
        </head>
        <body>
            {%app_entry%}
            <footer>
                {%config%}
                {%scripts%}
                {%renderer%}
            </footer>
        </body>
    </html>
    ''')

# Create initial figures with Inter font
fig_spatial = px.scatter(df_umap, x='x', y='y', 
                        color='Cell Class',
                        color_discrete_sequence=px.colors.qualitative.D3)
fig_spatial.update_layout(
    dragmode='select', 
    title=dict(
        text='<b>Spatial Coordinates</b>',
        x=0.5,
        xanchor='center',
        font=dict(size=20, family='Inter')
    ),
    height=600,
    font=dict(family='Inter'),
    template='simple_white',
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)'
)

fig_umap = px.scatter(df_umap, x='UMAP1', y='UMAP2', 
                      color='Cell Class',
                      color_discrete_sequence=px.colors.qualitative.D3)
fig_umap.update_layout(
    title=dict(
        text='<b>UMAP Projection</b>',
        x=0.5,
        xanchor='center',
        font=dict(size=20, family='Inter')
    ),
    height=400,
    font=dict(family='Inter'),
    template='simple_white',
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)'
)

# Create initial bubble plot
initial_bubble_df = bubble_df_dict[('UMAP', n_clusters_options[0])]
fig_bubble = px.scatter(initial_bubble_df, x='Cluster', y='Gene', 
                       size='Mean Expression', color='Mean Expression',
                       size_max=15, title='Top Genes per Cluster')
fig_bubble.update_layout(
    title=dict(
        text='<b>Top Genes per Cluster</b>',
        x=0.5,
        xanchor='center',
        font=dict(size=20, family='Inter')
    ),
    height=800, 
    margin=dict(l=50, r=50, t=70, b=50),
    font=dict(family='Inter'),
    template='simple_white',
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)',
    coloraxis_colorbar=dict(
        thickness=10,
        len=0.5,
        title=dict(
            text='Mean<br>Expression',
            font=dict(family='Inter')
        )
    )
)

# Define app layout with Inter font
app.layout = html.Div([
    html.H1("♠ SPADE", 
            style={'textAlign': 'left', 'margin-bottom': '10px', 'fontFamily': 'Inter', 'fontSize': '24px', 'padding-left': '20px'}),
    
    html.Div([
        # Controls column
        html.Div([
            html.Div([
                html.P("SPADE (Simultaneous Projection and Analysis of Dimensionality-reduced Expression) is an interactive visualization tool designed for exploring and analyzing high-dimensional spatial transcriptomic data.",
                       style={'textAlign': 'left', 'margin-bottom': '20px', 'fontFamily': 'Inter', 'fontSize': '14px', 'line-height': '1.5'}),

                
                html.Div([
                    html.Label("Dimension Reduction Type", style={'font-weight': 'bold', 'fontFamily': 'Inter', 'fontSize': '16px'}),
                    dcc.RadioItems(
                        id='dim-reduction-method',
                        options=[
                            {'label': 'UMAP', 'value': 'umap'},
                            {'label': 't-SNE', 'value': 'tsne'}
                        ],
                        value='umap',
                        inline=True,
                        style={'margin-top': '10px', 'margin-bottom': '20px', 'fontFamily': 'Inter', 'fontSize': '14px'}
                    ),
                ]),
                
                html.Div([
                    html.Label("Number of Clusters", style={'font-weight': 'bold', 'fontFamily': 'Inter', 'fontSize': '16px'}),
                    html.Div(
                        dcc.Slider(
                            id='cluster-number',
                            min=min(n_clusters_options),
                            max=max(n_clusters_options),
                            step=1,
                            marks={i: str(i) for i in n_clusters_options},
                            value=n_clusters_options[0],
                            tooltip={'placement': 'bottom'}
                        ),
                        style={'margin-top': '10px', 'fontFamily': 'Inter', 'fontSize': '14px'}
                    ),
                ], style={'margin-top': '20px'}),
                
                html.Div([
                    html.Label("Selection Controls", style={'font-weight': 'bold', 'fontFamily': 'Inter', 'fontSize': '16px', 'margin-top': '20px', 'display': 'block'}),
                    html.Button(
                        'Reset Selection',
                        id='reset-button',
                        n_clicks=0,
                        style={
                            'margin-top': '10px',
                            'padding': '10px 20px',
                            'background-color': '#f8f9fa',
                            'border': '1px solid #ddd',
                            'border-radius': '5px',
                            'font-family': 'Inter',
                            'font-weight': 'bold',
                            'fontSize': '14px',
                            'cursor': 'pointer',
                            'display': 'block'
                        }
                    )
                ])
            ], style={'padding': '20px', 'background-color': '#f5f5f5', 'border-radius': '5px'})
        ], style={'width': '15%', 'display': 'inline-block', 'vertical-align': 'top', 
                 'padding-right': '20px', 'background-color': '#f5f5f5', 'border-radius': '5px'}),
        
        # Visualizations column
        html.Div([
            dcc.Graph(id='spatial-graph', figure=fig_spatial),
            dcc.Graph(id='umap-graph', figure=fig_umap)
        ], style={'width': '40%', 'display': 'inline-block', 'vertical-align': 'top',
                 'padding': '10px'}),
        # Bubble plot and cluster visualization column
        html.Div([
            html.Div([
                html.Div([
                    dcc.Graph(id='bubble-plot', figure=fig_bubble, style={'height': '600px'})
                ], style={'width': '60%', 'display': 'inline-block', 'vertical-align': 'top'}),
                
                html.Div([
                    html.Label("Select Cluster for Cell Type Composition", style={'font-weight': 'bold', 'fontFamily': 'Inter', 'fontSize': '16px'}),
                    dcc.Dropdown(
                        id='cluster-selector',
                        options=[{'label': f'Cluster {i}', 'value': str(i)} for i in range(n_clusters_options[0])],
                        value='0',
                        style={'width': '100%', 'margin-bottom': '10px', 'fontFamily': 'Inter'}
                    ),
                    dcc.Graph(id='pie-chart', figure=fig_pie, style={'height': '500px'})
                ], style={'width': '40%', 'display': 'inline-block', 'vertical-align': 'top', 'padding-left': '20px'})
            ], style={'display': 'flex', 'margin-bottom': '20px', 'height': '650px'}),
            
            html.Div([
                dcc.Graph(id='cluster-graph', figure={}, style={'height': '400px'})
            ], style={'margin-bottom': '20px'})
        ], style={'width': '45%', 'display': 'inline-block', 'vertical-align': 'top',
                 'padding': '10px', 'border-left': '1px solid #ddd'})
    ], style={'display': 'flex', 'border': '1px solid #ddd', 'border-radius': '5px'})
])

# Callbacks
@app.callback(
    Output('cluster-graph', 'figure'),
    Input('dim-reduction-method', 'value'),
    Input('cluster-number', 'value')
)
def update_cluster_graph(dim_method, n_clusters):
    if dim_method == 'umap':
        x_col, y_col = 'UMAP1', 'UMAP2'
        title_prefix = 'UMAP'
        df = df_umap
        cluster_col = f'umap_cluster_{n_clusters}'
    else:
        x_col, y_col = 'tSNE1', 'tSNE2'
        title_prefix = 't-SNE'
        df = df_tsne
        cluster_col = f'tsne_cluster_{n_clusters}'
    
    fig = px.scatter(
        df, x=x_col, y=y_col, 
        color=cluster_col,
        color_discrete_sequence=px.colors.qualitative.D3,
        title=f"{title_prefix} with {n_clusters} Clusters",
        labels={x_col: f'{title_prefix}1', y_col: f'{title_prefix}2'},
        category_orders={cluster_col: sorted(df[cluster_col].unique())},
        hover_data=['Cell Class']
    )
    
    fig.update_layout(
        title=dict(
            text=f'<b>{title_prefix} with {n_clusters} Clusters</b>',
            x=0.5,
            xanchor='center',
            font=dict(size=20, family='Inter')
        ),
        dragmode='select',
        legend_title='Cluster',
        height=400,
        font=dict(family='Inter'),
        template='simple_white',
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    
    return fig

@app.callback(
    Output('spatial-graph', 'figure'),
    Input('reset-button', 'n_clicks'),
    prevent_initial_call=True
)
def reset_selection(n_clicks):
    if n_clicks > 0:
        fig_spatial.update_traces(
            selectedpoints=None,
            unselected=dict(marker=dict(opacity=1))
        )
    return fig_spatial

@app.callback(
    Output('umap-graph', 'figure'),
    [Input('dim-reduction-method', 'value'),
     Input('spatial-graph', 'selectedData'),
     Input('reset-button', 'n_clicks')],
    prevent_initial_call=True
)
def update_umap_graph(dim_method, selectedData, n_clicks):
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update
        
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if dim_method == 'umap':
        df = df_umap
        x_col, y_col = 'UMAP1', 'UMAP2'
        title = 'UMAP Projection'
    else:
        df = df_tsne
        x_col, y_col = 'tSNE1', 'tSNE2'
        title = 't-SNE Projection'
    
    # Create a new figure instance
    fig = px.scatter(df, x=x_col, y=y_col, 
                    color='Cell Class',
                    color_discrete_sequence=px.colors.qualitative.D3)
    
    fig.update_layout(
        title=dict(
            text=f'<b>{title}</b>',
            x=0.5,
            xanchor='center',
            font=dict(size=20, family='Inter')
        ),
        height=400,
        font=dict(family='Inter'),
        template='simple_white',
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    
    if trigger_id == 'reset-button' and n_clicks:
        # Clear selection and reset opacity
        fig.update_traces(
            mode='markers',
            selectedpoints=None,
            opacity=1
        )
    elif selectedData and 'points' in selectedData:
        selected_indices = [point['pointIndex'] for point in selectedData['points']]
        # Update opacity for selected and unselected points
        fig.update_traces(
            mode='markers',
            selectedpoints=selected_indices,
            selected=dict(marker=dict(color='red', size=10)),
            unselected=dict(marker=dict(opacity=0.3))
        )
    
    return fig

@app.callback(
    Output('bubble-plot', 'figure'),
    Input('dim-reduction-method', 'value'),
    Input('cluster-number', 'value')
)
def update_bubble_plot(dim_method, n_clusters):
    bubble_df = bubble_df_dict[(dim_method.upper(), n_clusters)]
    
    fig = px.scatter(bubble_df, x='Cluster', y='Gene', 
                    size='Mean Expression', color='Mean Expression',
                    size_max=15, title='Top Genes per Cluster')
    
    fig.update_layout(
        title=dict(
            text='<b>Top Genes per Cluster</b>',
            x=0.5,
            xanchor='center',
            font=dict(size=20, family='Inter')
        ),
        height=600,
        margin=dict(l=50, r=50, t=70, b=50),
        xaxis_title='Cluster',
        yaxis_title='Gene',
        showlegend=False,
        font=dict(family='Inter'),
        template='simple_white',
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        coloraxis_colorbar=dict(
            thickness=10,
            len=0.5,
            title=dict(
                text='Mean<br>Expression',
                font=dict(family='Inter')
            )
        ),
        yaxis=dict(
            automargin=True,
            tickfont=dict(size=10)
        ),
        xaxis=dict(
            automargin=True,
            tickfont=dict(size=10)
        )
    )
    
    return fig

@app.callback(
    Output('cluster-selector', 'options'),
    Output('cluster-selector', 'value'),
    Input('dim-reduction-method', 'value'),
    Input('cluster-number', 'value')
)
def update_cluster_selector(dim_method, n_clusters):
    # Update dropdown options based on number of clusters
    options = [{'label': f'Cluster {i}', 'value': str(i)} for i in range(n_clusters)]
    # Reset to first cluster
    value = '0'
    return options, value

@app.callback(
    Output('pie-chart', 'figure'),
    Input('dim-reduction-method', 'value'),
    Input('cluster-number', 'value'),
    Input('cluster-selector', 'value')
)
def update_pie_chart(dim_method, n_clusters, selected_cluster):
    method_key = 'UMAP' if dim_method == 'umap' else 'TSNE'
    cell_counts = cell_composition_dict[(method_key, n_clusters)][selected_cluster]
    
    return create_pie_chart(
        cell_counts, 
        f'Cell Type Composition<br>{method_key} Cluster {selected_cluster}'
    )

if __name__ == '__main__':
    app.run(debug=True) 