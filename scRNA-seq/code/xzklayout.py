#%%
import os
os.chdir('F:\\LabW')
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import pandas as pd
import numpy as np 
import dash
from xzkapp import app
#%%
# formatting info
corporate_colors = {
    'dark-blue-grey' : 'rgb(62, 64, 76)',
    'medium-blue-grey' : 'rgb(77, 79, 91)',
    'superdark-green' : 'rgb(41, 56, 55)',
    'dark-green' : 'rgb(57, 81, 85)', #################
    'medium-green' : 'rgb(93, 113, 120)',
    'light-green' : 'rgb(186, 218, 212)',
    'pink-red' : 'rgb(255, 101, 131)',
    'dark-pink-red' : 'rgb(247, 80, 99)',
    'white' : 'rgb(251, 251, 252)',
    'light-grey' : 'rgb(208, 206, 206)',
    'azure3':'RGB(193,205,205)',
    'sgigray92':'RGB(183,183,183)',
    'whitesmoke':'RGB(245,245,245)'
}

externalgraph_rowstyling = {
    # 'margin-left' : '15px',
    # 'margin-right' : '15px',
    'background-color' : corporate_colors['whitesmoke'],
}

externalgraph_colstyling = {
    'border-radius' : '10px',
    'border-style' : 'solid',
    'border-width' : '1px',
    'border-color' : corporate_colors['azure3'],
    'background-color' : corporate_colors['azure3'],
    'box-shadow' : '0px 0px 17px 0px rgba(186, 218, 212, .5)',
    'padding-top' : '10px'
}

dropdown1div_borderstyling = {
    'border-radius' : '0px 0px 0px 0px',
    'border-style' : 'solid',
    'border-width' : '1px',
    'border-color' : corporate_colors['light-green'],
    'background-color' : corporate_colors['light-green'],
    # 'box-shadow' : '2px 5px 5px 1px rgba(255, 101, 131, .5)'
    }

dropdown2div_borderstyling = {
    'border-radius' : '0px 0px 0px 0px',
    'border-style' : 'solid',
    'border-width' : '1px',
    'border-color' : corporate_colors['light-green'],
    'background-color' : corporate_colors['light-green'],
    # 'box-shadow' : '2px 5px 5px 1px rgba(255, 101, 131, .5)'
    }

dropdown3div_borderstyling = {
    'border-radius' : '0px 0px 10px 10px',
    'border-style' : 'solid',
    'border-width' : '1px',
    'border-color' : corporate_colors['light-green'],
    'background-color' : corporate_colors['light-green'],
    'box-shadow' : '2px 5px 5px 1px rgba(255, 101, 131, .5)'
    }

# import tnse data
UCLA1_tsne_data = pd.read_csv('./20210323tsne_data1.csv', index_col=0)
UCLA2_tsne_data = pd.read_csv('./20210323tsne_data2.csv', index_col=0)
UCLA1_tsne_data['cluster'] = UCLA1_tsne_data['cluster'].map(lambda x: "cluster_" + str(x))
UCLA2_tsne_data['cluster'] = UCLA2_tsne_data['cluster'].map(lambda x: "cluster_" + str(x))

#tsne_data = pd.read_csv('./20210116tsne_dataGG.csv',index_col=0)
#print(tsne_data["cluster"])
#tsne_data["cluster"] = tsne_data['cluster'].map(lambda x: "cluster_" + str(x))
#print(tsne_data["cluster"])

# import gene name data
UCLA1_genes_data = pd.read_csv('./xzk20210323rna_data1.csv')
UCLA2_genes_data = pd.read_csv('./xzk20210323rna_data2.csv')
#%%
# import human genes correlated with NANOS3 data
UCLA1_genes_NANOS3 = pd.read_csv('./UCLA1_NANOS3_cor_allgene.csv')
UCLA2_genes_NANOS3 = pd.read_csv('./UCLA2_NANOS3_cor_allgene.csv')
UCLA1_genes_NANOS3.head()
UCLA1_genes_NANOS3.iloc[:,0]
#%%

# create dataset dropdown
dataset_dropdown = [{"label": 'UCLA1', "value": 'UCLA1'}, 
                    {"label": 'UCLA2', "value": 'UCLA2'}]

# create gene dropdown
UCLA1_genes = UCLA1_genes_NANOS3.iloc[:,0].to_list()
#UCLA1_genes_dropdown = [
#    {"label":i, "value":i} for i in UCLA1_genes
# ]

UCLA2_genes = UCLA2_genes_NANOS3.iloc[:,0].to_list()
#UCLA2_genes_dropdown = [
#    {"label":i, "value":i} for i in UCLA2_genes
#]
all_options = {
    'UCLA1': UCLA1_genes,
    'UCLA2': UCLA2_genes
}

# Header with logo
def get_header():
    header = html.Div([
        html.Div([], className = 'col-4'),
        html.Div([
                html.H1(children = "Human PGCLC Specific Gene", 
                style = {'textAligh':'center','color':'white'})
                ],
                className = 'col-5',
                style = {'padding-top':'1%'}
        ),
        html.Div([html.Img(src = app.get_asset_url('labw.png'),
                height = '43 px',
                width = 'auto')],
                    className = "col-1",
                    style = {'align-items': 'center',
                    'padding-top' : '1%',
                    'height' : 'auto'})
        ], 
        className = 'row',
        style = {'height': '3%', 'background-color': corporate_colors['superdark-green']}
        )
    return(header)

# get footer
def get_footer():
    footer = html.Div([
        html.Div([], className='col-3'),
        html.Div([
            html.Br(),
            html.Br(),
            html.P(children = 'Mailing address', 
                    style = {'textAlign' : 'center','margin-bottom':0}
                    ),
            html.P(children = 'A314, ZJU-UoE Institute,',
                    style = {'textAlign' : 'center','margin-bottom':0}
                    ),
            html.P(children = '718 East Haizhou Rd.,',
                    style = {'textAlign' : 'center','margin-bottom':0}
                    ),
            html.P(children = 'Haining, Zhejiang 314400,',
                    style = {'textAlign' : 'center','margin-bottom':0}
                    ),
            html.P(children = 'P.R. China',
                    style = {'textAlign' : 'center','margin-bottom':0}
                    ),        
            html.Br(),
            html.Br(),
            html.Br()
            ], className='col-6'),
        html.Div([], className = 'col-3')      
    ], className = 'row',    
        style = {'height': '5%', 'background-color': corporate_colors['whitesmoke']}
    )
    return(footer)

# empty row
def get_emptyrow(**mystyle):
    emptyrow = html.Div([
            html.Div([
            html.Br()
        ], className = 'col-12')
    ],
    className = 'row',
    style = mystyle)
    return(emptyrow)


# graph
totals = html.Div([
    get_header(), # row 1
    get_emptyrow(**{'height':'30px', #row 2
                    'background-color' : corporate_colors['light-green']}), # Row 2
    html.Div([ # external row
        html.Div([], className= 'col-1'), # blank internal columns
        html.Div([ 
            html.Div([# Dataset dropdown
                html.H2(children= 'Dataset'), 
                dcc.Dropdown(
                    id= 'dataset_dropdown',
                    options= [{'label': k, 'value': k} for k in all_options.keys()],
                    value= 'UCLA1',
                    multi= False,
                    style= {'font-size': '13px', 
                            'color' : corporate_colors['medium-blue-grey'], 
                            'white-space': 'nowrap', 
                            'text-overflow': 'ellipsis'}
                )], className= 'col-6', style = {'float': 'left'}),
            html.Div([# human genes correlated with NANOS3 dropdown
                html.H2(children= 'Human Genes'),
                dcc.Dropdown(
                    id= 'human_gene_correlated_with_NANOS3',
                    #options= [{'label': j, 'value': j} for j in all_options.values()], ###########??????????
                    value= 'NANOS3',
                    multi= False,
                    style= {'font-size': '13px', 'color' : corporate_colors['medium-blue-grey'], 'white-space': 'nowrap', 'text-overflow': 'ellipsis'}
                )], className= 'col-6', style = {'float': 'right'}),
            ]),
        html.Div([
            html.Br(),
            html.Div(children='Correlation: ', id='correlation', className= 'col-6', style = {'float': 'right'})
            ], className = 'row'),
    ]),

    # empty row 3
    get_emptyrow(**{'height':'30px','background-color' : corporate_colors['light-green']}),
    
    # Row4 graphs
    html.Div([ # External row
        html.Div([ # Internal row
            html.Div([], className= 'col-1'),# blank 1 row
            html.Div([
                dcc.Graph(id='UMAP')
            ], className='col-5', style = {'float': 'left'}), 
            html.Div([
                dcc.Graph(id='NANOS3_graph')
            ], className='col-5', style = {'float': 'right'}),
            
        html.Div([], className='col-1') # Blank 1 column
        ], className= 'col-10')
    ], className = 'row'), 
    # row5 graphs
    html.Div([ # External row
        html.Div([ # Internal row
            html.Div([], className= 'col-1'),# blank 1 row
            html.Div([
                dcc.Graph(id='expression_level_graph')
            ], className='col-5', style = {'float': 'left'}),
            html.Div([
                dcc.Graph(id='curve_plot_graph')
            ], className='col-5', style = {'float': 'right'})
        ], className= 'col-10')
    ], className = 'row'), 
    get_footer()
])


