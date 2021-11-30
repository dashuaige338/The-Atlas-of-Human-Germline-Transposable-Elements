#%%
import os
os.chdir('F:\\LabW')
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
import pandas as pd
#%%
import numpy as np
from xzkapp import app
from xzklayout import all_options
from dash.dependencies import Input, Output, State
import plotly.express as px
import dash
from scipy import sparse
#%%
####################### Corporate css formatting
corporate_colors = {
    'dark-blue-grey' : 'rgb(62, 64, 76)',
    'medium-blue-grey' : 'rgb(77, 79, 91)',
    'superdark-green' : 'rgb(41, 56, 55)',
    'dark-green' : 'rgb(57, 81, 85)',
    'medium-green' : 'rgb(93, 113, 120)',
    'light-green' : 'rgb(186, 218, 212)',
    'pink-red' : 'rgb(255, 101, 131)',
    'dark-pink-red' : 'rgb(247, 80, 99)',
    'white' : 'rgb(251, 251, 252)',
    'light-grey' : 'rgb(208, 206, 206)',
    'ghostwhite':'RGB(248,248,255)'
}

####################### Corporate chart formatting

corporate_title = {
    'font' : {
        'size' : 16,
        'color' : corporate_colors['white']}
}

corporate_xaxis = {
    'showgrid' : False,
    'linecolor' : corporate_colors['dark-green'],
    'color' : corporate_colors['dark-green'],
    'tickangle' : 315,
    'titlefont' : {
        'size' : 12,
        'color' : corporate_colors['dark-green']},
    'tickfont' : {
        'size' : 11,
        'color' : corporate_colors['dark-green']},
    'zeroline': False
}

corporate_yaxis = {
    'showgrid' : False,
    'color' : corporate_colors['dark-green'],
    # 'gridwidth' : 0.5,
    # 'gridcolor' : corporate_colors['dark-green'],
    'linecolor' : corporate_colors['dark-green'],
    'titlefont' : {
        'size' : 12,
        'color' : corporate_colors['dark-green']},
    'tickfont' : {
        'size' : 11,
        'color' : corporate_colors['dark-green']},
    'zeroline': False
}

corporate_font_family = 'Arial'

corporate_legend = {
    'orientation' : 'v',
    'yanchor' : 'top',
    'y' : 1.21,
    'xanchor' : 'right',
    'x' : 1.05,
	'font' : {'size' : 12, 'color' : corporate_colors['light-grey']}
} # Legend will be on the top right, above the graph, horizontally

corporate_margins = {'l' : 5, 'r' : 5, 't' : 45, 'b' : 15}  # Set top margin to in case there is a legend

corporate_layout = go.Layout(
    font = {'family' : corporate_font_family},
    title = corporate_title,
    title_x = 0.5, # Align chart title to center
    paper_bgcolor = 'rgba(0,0,0,0)',
    plot_bgcolor = 'rgba(0,0,0,0)',
    xaxis = corporate_xaxis,
    yaxis = corporate_yaxis,
    # height = 270,
    # legend = corporate_legend,
    margin = corporate_margins
)



############ import data ############
## import tsne data
UCLA1_tsne_data = pd.read_csv('./20210323tsne_data1.csv', index_col=0)
UCLA2_tsne_data = pd.read_csv('./20210323tsne_data2.csv', index_col=0)
UCLA1_tsne_data['cluster'] = UCLA1_tsne_data['cluster'].map(lambda x: "cluster_" + str(x))
UCLA2_tsne_data['cluster'] = UCLA2_tsne_data['cluster'].map(lambda x: "cluster_" + str(x))
UCLA1_data_cluster = UCLA1_tsne_data['cluster'].to_list()
UCLA2_data_cluster = UCLA2_tsne_data['cluster'].to_list()
UCLA1_data_cluster = list(sorted(set(UCLA1_data_cluster), key=lambda x: int(x.split('_')[1])))
UCLA2_data_cluster = list(sorted(set(UCLA2_data_cluster), key=lambda x: int(x.split('_')[1])))
UCLA1_data_color = px.colors.qualitative.Dark24[:len(UCLA1_data_cluster)]
UCLA2_data_color = px.colors.qualitative.Dark24[:len(UCLA2_data_cluster)]
#%%
## import gene data
UCLA1_gene_expression_file = './xzk20210323rna_data1.csv'
UCLA2_gene_expression_file = './xzk20210323rna_data2.csv'
UCLA1_gene_data = pd.read_csv(UCLA1_gene_expression_file)
UCLA2_gene_data = pd.read_csv(UCLA2_gene_expression_file)
#%%
## import correlation-NANOS3 data
UCLA1_genes_NANOS3 = pd.read_csv('./UCLA1_NANOS3_cor_allgene.csv')
UCLA2_genes_NANOS3 = pd.read_csv('./UCLA2_NANOS3_cor_allgene.csv')
UCLA1_genes_NANOS3.rename(columns={'Unnamed: 0': 'gene', 'V1':'correlation'}, inplace=True)
UCLA2_genes_NANOS3.rename(columns={'Unnamed: 0': 'gene', 'V1':'correlation'}, inplace=True)
UCLA1_genes_NANOS3.set_index(UCLA1_genes_NANOS3['gene'], inplace=True)
UCLA2_genes_NANOS3.set_index(UCLA2_genes_NANOS3['gene'], inplace=True)
############################## plot #############################
# plot tSNE function
def plottSNE(group, data, tempColormap, plottype):
    cluster = {
        'UCLA1': UCLA1_data_cluster,
        'UCLA2': UCLA2_data_cluster
    }
    fig = px.scatter(data, x = 'UMAP_1',y = 'UMAP_2', color= plottype,
                color_discrete_map=tempColormap,
                category_orders={"cluster": cluster[group]},
                width=400, height=400
                )
    return fig

# featureplot function
def plotFeaturePlot(group, gene):
    expression_file_route = {
        'UCLA1': UCLA1_gene_expression_file,
        'UCLA2': UCLA2_gene_expression_file
    }
    all_tsne_data = {
        'UCLA1': UCLA1_tsne_data, 
        'UCLA2': UCLA2_tsne_data
    }
    expression_data = pd.read_csv(expression_file_route[group], usecols= [gene])[gene].to_list()
    all_tsne_data[group]['expression_level'] = expression_data
    featureplot = px.scatter(all_tsne_data[group], x='UMAP_1', y='UMAP_2',
            color='expression_level',
            color_continuous_scale=px.colors.sequential.BuPu,
            width=400, height=400
        )
    return featureplot
#%%
# curve plot
# expression level
def expression_mean(type, dataset):
    if dataset == 'UCLA1':
        dataindex = list(UCLA1_gene_data.iloc[:,0])
        a = []
        for i in dataindex:
            a.append(type in i)
        cell_mean = UCLA1_gene_data.loc[a].mean(axis=0).to_frame(name=type)
        cell_error = UCLA1_gene_data.loc[a].sem(axis=0).to_frame(name='error_'+type)
        cell = cell_mean.join(cell_error)
    elif dataset == 'UCLA2':
        dataindex = list(UCLA2_gene_data.iloc[:,0])
        a = []
        for i in dataindex:
            a.append(type in i)
        cell_mean = UCLA2_gene_data.loc[a].mean(axis=0).to_frame(name=type)
        cell_error = UCLA2_gene_data.loc[a].sem(axis=0).to_frame(name='error_'+type)
    return cell_mean, cell_error

UCLA1_U1hESC = expression_mean('hESC', 'UCLA1')
UCLA1_U1iMeLC = expression_mean('iMeLC', 'UCLA1')
UCLA1_U1d1 = expression_mean('d1', 'UCLA1')
UCLA1_U1d2 = expression_mean('d2', 'UCLA1')
UCLA1_U1d3 = expression_mean('d3', 'UCLA1')
UCLA1_U1d4 = expression_mean('d4', 'UCLA1')
UCLA1_expression_level = UCLA1_U1hESC[0].join([UCLA1_U1iMeLC[0], UCLA1_U1d1[0], UCLA1_U1d2[0], UCLA1_U1d3[0], UCLA1_U1d4[0]])
UCLA1_expression_level_error = UCLA1_U1hESC[1].join([UCLA1_U1iMeLC[1], UCLA1_U1d1[1], UCLA1_U1d2[1], UCLA1_U1d3[1], UCLA1_U1d4[1]])

UCLA2_U2hESC = expression_mean('hESC', 'UCLA2')
UCLA2_U2iMeLC = expression_mean('iMeLC', 'UCLA2')
UCLA2_U2d1 = expression_mean('d1', 'UCLA2')
UCLA2_U2d2 = expression_mean('d2', 'UCLA2')
UCLA2_U2d3 = expression_mean('d3', 'UCLA2')
UCLA2_U2d4 = expression_mean('d4', 'UCLA2')
UCLA2_expression_level = UCLA2_U2hESC[0].join([UCLA2_U2iMeLC[0], UCLA2_U2d1[0], UCLA2_U2d2[0], UCLA2_U2d3[0], UCLA2_U2d4[0]])
UCLA2_expression_level_error = UCLA2_U2hESC[1].join([UCLA2_U2iMeLC[1], UCLA2_U2d1[1], UCLA2_U2d2[1], UCLA2_U2d3[1], UCLA2_U2d4[1]])

#%%
def curve_plot(dataset, gene):
    meangroup = {
        'UCLA1': UCLA1_expression_level,
        'UCLA2': UCLA2_expression_level
    }
    errorgroup = {
        'UCLA1': UCLA1_expression_level_error,
        'UCLA2': UCLA2_expression_level_error
    }
    x = ['hESC', 'hiMeLC', 'D1', 'D2', 'D3', 'D4']
    fig = go.Figure(data=go.Scatter(
        y=list(meangroup[dataset].loc[gene,:]),
        x=x,
        error_y=dict(type='data', array=list(errorgroup[dataset].loc[gene,:]), visible=True)
        ),
        layout=go.Layout(height=400, width=400)
    )
    return fig

def correlation(dataset, gene):
    corgroup = {
        'UCLA1': UCLA1_genes_NANOS3,
        'UCLA2': UCLA2_genes_NANOS3
    }
    cor = str(corgroup[dataset].loc[gene, 'correlation'])
    return cor

######################### callback #########################
#%%
####### dropdown callback
@app.callback(
    Output("human_gene_correlated_with_NANOS3", "options"),
    Input("dataset_dropdown", 'value')
)
def callback_dropdown(selected_dataset):
    return [{'label': i, 'value': i} for i in all_options[selected_dataset]]

####### tsne callback
@app.callback(
    dash.dependencies.Output(component_id='UMAP', component_property= 'figure'),
    [dash.dependencies.Input('dataset_dropdown', 'value')]
)
def update_cluster_graph(dataset):
    UCLA1_data_cluster_color = dict(zip(UCLA1_data_cluster, UCLA1_data_color))
    UCLA2_data_cluster_color = dict(zip(UCLA2_data_cluster, UCLA2_data_color))
    data_cluster_color = {
        'UCLA1':UCLA1_data_cluster_color, 
        'UCLA2':UCLA2_data_cluster_color
    }
    tsne_data = {
        'UCLA1': UCLA1_tsne_data, 
        'UCLA2': UCLA2_tsne_data
    }
    fig = plottSNE(dataset, tsne_data[dataset], data_cluster_color[dataset], 'cluster')

    fig.update_layout(corporate_layout)
    return(fig)

###### featureplot callback
# NANOS3 featureplot
@app.callback(
    dash.dependencies.Output('NANOS3_graph', 'figure'), 
    dash.dependencies.Input('dataset_dropdown', 'value')
)
def update_human_NANOS3_featureplot(dataset):
    fig = plotFeaturePlot(dataset, 'NANOS3')
    fig.update_layout(corporate_layout)
    return fig

# featureplot
@app.callback(
    dash.dependencies.Output('expression_level_graph', 'figure'),
    dash.dependencies.Input('human_gene_correlated_with_NANOS3', 'value'), 
    dash.dependencies.Input('dataset_dropdown', 'value')
)
def update_human_gene_featureplot(gene='human_gene_correlated_with_NANOS3', dataset='dataset_dropdown'):
    fig = plotFeaturePlot(dataset, gene)
    fig.update_layout(corporate_layout)
    return fig

# curve plot
@app.callback(
    dash.dependencies.Output('curve_plot_graph', 'figure'),
    dash.dependencies.Input('human_gene_correlated_with_NANOS3', 'value'), 
    dash.dependencies.Input('dataset_dropdown', 'value')
)
def update_human_gene_curveplot(gene='human_gene_correlated_with_NANOS3', dataset='dataset_dropdown'):
    fig = curve_plot(dataset, gene)
    fig.update_layout(corporate_layout)
    return fig

# correlation
@app.callback(
    dash.dependencies.Output('correlation', 'children'),
    dash.dependencies.Input('human_gene_correlated_with_NANOS3', 'value'), 
    dash.dependencies.Input('dataset_dropdown', 'value')
)
def update_correlation(gene='human_gene_correlated_with_NANOS3', dataset='dataset_dropdown'):
    corr = correlation(dataset, gene)
    return "Correlation with NANOS3: {}".format(corr)