from flask import Flask, render_template, request, flash, redirect, url_for, session
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.serving import run_simple
import base64
import datetime
import io
import numpy as np
import pandas as pd
import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import plotly.express as px
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import os
import psycopg2
import json

flaskApp = Flask(__name__)
flaskApp.secret_key = os.urandom(24)

@flaskApp.route('/hPGCLCdb', methods=["GET"])
def tsne():
    return render_template("index.html")

############
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app1 = dash.Dash(__name__, server=flaskApp, external_stylesheets=[dbc.themes.DARKLY],routes_pathname_prefix='/hPGCLCdb/bulk/', requests_pathname_prefix='/hPGCLCdb/bulk/')
app1.title = "hPGCLCdb"

TEdata1 = pd.read_csv('./bulk_RNA_seq/20210311Xinyu_TEexpression_RPKM_onlyTEwithReads_1180_detailedinfo.csv')
TE_id = TEdata1['type'].tolist()

TEdata2 = TEdata1.drop(['type','Length'],axis = 1)

TE_mean=pd.DataFrame(TEdata1,columns=['type','hESC_mean','iMeLC_mean','PGC_mean','PGCLC_mean'])

sample = list(TEdata2)

genedata = pd.read_csv('./bulk_RNA_seq/gene.csv')
symbol = genedata['symbol']
gene_mean=pd.DataFrame(genedata,columns=['symbol','hESC_mean','iMeLC_mean','PGC_mean','PGCLC_mean'])



app1.layout = html.Div([
    html.Div([
    html.Br(),
    html.H3('Samples',style={'textAlign': 'center'}),
    html.Div([dcc.Dropdown(
                  id = 'sample_name',
                  options = [{'label':i, 'value':i} for i in sample],
                  value = '',multi=True,style={"color":"black"})]),
    html.Br(),
    html.H3('Genes',style={'textAlign': 'center'}),
    html.Div([dcc.Dropdown(
                  id = 'gene_name',
                  options = [{'label':i, 'value':i} for i in symbol],
                  value = '',multi=True,style={"color":"black"})]),
    html.Br(),
    html.H3('Transposable Elements (TE)',style={'textAlign': 'center'}),
    html.Div([dcc.Dropdown(
                  id = 'TE_name',
                  options = [{'label':i, 'value':i} for i in TE_id],
                  value = '',multi=True,style={"color":"black"})]),
    ],className='col'),
    html.Div([
      html.Div([
        html.Br(),
        html.Br(),
        html.H3('Heatmap for gene expression levels',style={'textAlign': 'center'}),
        html.Br(),
        html.Div([dcc.Graph(id='heatmap1')]),
        html.Br(),
        html.Br(),
        html.H3('gene expression change with stages',style={'textAlign': 'center'}),
        html.Br(),
        html.Div([dcc.Graph(id='linechart1')]),
      ],className='col'),
      html.Div([
        html.Br(),
        html.Br(),
        html.H3('Heatmap for TE expression levels',style={'textAlign': 'center'}),
        html.Br(),
        html.Div([dcc.Graph(id='heatmap2')]),
        html.Br(),
        html.Br(),
        html.H3('TE expression change with stages',style={'textAlign': 'center'}),
        html.Br(),
        html.Div([dcc.Graph(id='linechart2')])
      ],className='col')    
],className='col')],className='col')


@app1.callback(Output('heatmap1','figure'),[Input('sample_name','value'),Input('gene_name','value')])
def update_figure1(sample_name,gene_name):
  TEtempDf = pd.DataFrame()
  if sample_name !='' and gene_name != '':
        for i in gene_name:
                  TEtempdata1 = genedata[(genedata['symbol'] == i)]
                  #print(TEtempdata1)
                  TEtempdata2 = TEtempdata1[sample_name]
                  #print(TEtempdata2)
                  TEtempDf = TEtempDf.append(pd.DataFrame(data=TEtempdata2))
  TEtempDf = TEtempDf.reset_index(drop=True)
  #print(TEtempDf)
  TEa = len(TEtempDf)
  #print(TEa)
  TEzdata1 = [TEtempDf.iloc[i] for i in range(TEa)]
  TEzdata = [np.log2(i+1) for i in TEzdata1]
  #print (TEzdata)
  TExdata = sample_name
  TEydata = gene_name

  # print('value:',value)
  trace = go.Heatmap(x=TExdata,y=TEydata,z=TEzdata,colorscale = 'Peach')
  fig = dict(
     data=[trace],
     layout = go.Layout(title= 'sample gene expression level (log2(RPKM+1))',
                       yaxis = dict (tickfont = dict(size=13)),
                       xaxis = dict(tickfont = dict (size=10)),
                       margin = dict(l = 110, b = 170, t = 60, r = 30)))
  return fig


@app1.callback(Output('linechart1','figure'),Input('gene_name','value'))
def update_figure2(gene_name):
  TEtempDf = pd.DataFrame()
  if gene_name != '':
        for i in gene_name:
                  TEtempdata1 = gene_mean[(gene_mean['symbol'] == i)]
                  #print(TEtempdata1)
                  #print(TEtempdata2)
                  TEtempDf = TEtempDf.append(pd.DataFrame(data=TEtempdata1))
  TEtempDf = TEtempDf.reset_index(drop=True)
  TEtempDf1 =pd.DataFrame(TEtempDf,columns=['hESC_mean','iMeLC_mean','PGC_mean','PGCLC_mean']) 
  #print(np.array(TEtempDf['type']).tolist())
  x=list(TEtempDf1)
  #print(x)
  TEa = len(TEtempDf1)
  data=[]
  #print(TEa)
  for i in range(TEa):
     y=np.array(TEtempDf1.iloc[i])
     y=np.array(y)
     tracei=go.Scatter(x=x,y=y,name=gene_name[i])
     data=data+[tracei]
  #print(data)
  #print (TEzdata)


  fig = dict(
     data=data,
     layout = go.Layout(title= 'changes of gene expression level with stage',
                       yaxis = dict (tickfont = dict(size=13)),
                       xaxis = dict(tickfont = dict (size=10)),
                       margin = dict(l = 110, b = 170, t = 60, r = 30)))
  return fig



@app1.callback(Output('heatmap2','figure'),[Input('sample_name','value'),Input('TE_name','value')])
def update_figure3(sample_name,TE_name):
  TEtempDf = pd.DataFrame()
  if sample_name !='' and TE_name != '':
        for i in TE_name:
                  TEtempdata1 = TEdata1[(TEdata1['type'] == i)]
                  #print(TEtempdata1)
                  TEtempdata2 = TEtempdata1[sample_name]
                  #print(TEtempdata2)
                  TEtempDf = TEtempDf.append(pd.DataFrame(data=TEtempdata2))
  TEtempDf = TEtempDf.reset_index(drop=True)
  #print(TEtempDf)
  TEa = len(TEtempDf)
  #print(TEa)
  TEzdata1 = [TEtempDf.iloc[i] for i in range(TEa)]
  TEzdata = [np.log2(i+1) for i in TEzdata1]
  #print (TEzdata)
  TExdata = sample_name
  TEydata = TE_name

  # print('value:',value)
  trace = go.Heatmap(x=TExdata,y=TEydata,z=TEzdata,colorscale = 'Peach')
  fig = dict(
     data=[trace],
     layout = go.Layout(title= 'sample TE expression level (log2(RPKM+1))',
                       yaxis = dict (tickfont = dict(size=13)),
                       xaxis = dict(tickfont = dict (size=10)),
                       margin = dict(l = 110, b = 170, t = 60, r = 30)))
  return fig


@app1.callback(Output('linechart2','figure'),Input('TE_name','value'))
def update_figure4(TE_name):
  TEtempDf = pd.DataFrame()
  if TE_name != '':
        for i in TE_name:
                  TEtempdata1 = TE_mean[(TE_mean['type'] == i)]
                  #print(TEtempdata1)
                  #print(TEtempdata2)
                  TEtempDf = TEtempDf.append(pd.DataFrame(data=TEtempdata1))
  TEtempDf = TEtempDf.reset_index(drop=True)
  TEtempDf1 =pd.DataFrame(TEtempDf,columns=['hESC_mean','iMeLC_mean','PGC_mean','PGCLC_mean']) 
  #print(np.array(TEtempDf['type']).tolist())
  x=list(TEtempDf1)
  #print(x)
  TEa = len(TEtempDf1)
  data=[]
  #print(TEa)
  for i in range(TEa):
     y=np.array(TEtempDf1.iloc[i])
     y=np.array(y)
     tracei=go.Scatter(x=x,y=y,name=TE_name[i])
     data=data+[tracei]
  #print(data)
  #print (TEzdata)


  fig = dict(
     data=data,
     layout = go.Layout(title= 'changes of TE expression level with stage',
                       yaxis = dict (tickfont = dict(size=13)),
                       xaxis = dict(tickfont = dict (size=10)),
                       margin = dict(l = 110, b = 170, t = 60, r = 30)))
  return fig

############

app2 = dash.Dash(__name__, server=flaskApp, external_stylesheets=[dbc.themes.DARKLY],routes_pathname_prefix='/hPGCLCdb/sc/', requests_pathname_prefix='/hPGCLCdb/sc/')
app2.title = "hPGCLCdb"

umap_data = pd.read_csv('./scRNA-seq/UCLA2_rep12.integrated.umap.csv')
#print(umap_data)
cluster1 = umap_data['Cluster']
cluster1 = np.array(cluster1).tolist()
cluster3 = list(set(cluster1))
#print(len(cluster3))
cluster2 = list(set(cluster1))
cluster2 = cluster2 + ['All']
#print(cluster2)
conn = psycopg2.connect("dbname = hpgcdb user = postgres password = Hrbl79eY host = 10.106.125.216 port = 5432")
c = conn.cursor()
c.execute("select name from expression")
gene0 = c.fetchall()
gene = []
for i in gene0:
   i = ''.join(i)
   gene = gene + [i]
#print(gene)
c.execute("select name from teexpression")
te0 = c.fetchall()
te = []
for i in te0:
   i = ''.join(i)
   te = te + [i]
#print(te)



app2.layout = html.Div([
  html.Div([
        html.Br(),
        html.H3('Choose clusters',style={'textAlign': 'center'}),
        html.Div([dcc.Dropdown(
                  id = 'cluster_name',
                  options = [{'label':i, 'value':i} for i in cluster2],
                  value = ['All'],multi=True,style={"color":"black"})]),
        html.Br(),
        html.H3('choose gene for featureplot',style={'textAlign': 'center'}),
        html.Div([dcc.Dropdown(
                  id = 'gene_name',
                  options = [{'label':i, 'value':i} for i in gene],
                  value = [''],multi=False,style={"color":"black"})]),
        html.Br(),
        html.H3('choose TE for featureplot',style={'textAlign': 'center'}),
        html.Div([dcc.Dropdown(
                  id = 'TE_name',
                  options = [{'label':i, 'value':i} for i in te],
                  value = [''],multi=False,style={"color":"black"})])
  ]),
  html.Div([
        html.Br(),
        html.Br(),
        html.H3('UMAP of gene expression and cluster',style={'textAlign': 'center'}),
        html.Br(),
        html.Div([dcc.Graph(id='umap')],style={'height':'600pt'}),
        html.Br(),
        html.H3('gene featureplot',style={'textAlign': 'center'}),
        html.Br(),
        html.Div([dcc.Graph(id='gene_featureplot')]),
        html.Br(),
        html.H3('TE featureplot',style={'textAlign': 'center'}),
        html.Br(),
        html.Div([dcc.Graph(id='TE_featureplot')]),
        html.Br(),
        html.H3('gene expression change with stages',style={'textAlign': 'center'}),
        html.Br(),
        html.Div([dcc.Graph(id='linechart1')]),
        html.Br(),
        html.H3('TE expression change with stages',style={'textAlign': 'center'}),
        html.Br(),
        html.Div([dcc.Graph(id='linechart2')])
  ])
],className = 'col')

def get_expression(gene_name):
    conn = psycopg2.connect("dbname = hpgcdb user = postgres password = Hrbl79eY host = 10.106.125.216 port = 5432")
    c = conn.cursor()
    c.execute("select expression from expression where name = '{}'".format(gene_name))
    results=c.fetchall()
    #print(results)
    result = results[0][0].tobytes()
    result = json.loads(result)
    result = pd.DataFrame(result,columns=['expression'])
    return result

def get_teexpression(gene_name):
    conn = psycopg2.connect("dbname = hpgcdb user = postgres password = Hrbl79eY host = 10.106.125.216 port = 5432")
    c = conn.cursor()
    c.execute("select expression from teexpression where name = '{}'".format(gene_name))
    results=c.fetchall()
    #print(results)
    result = results[0][0].tobytes()
    result = json.loads(result)
    result = pd.DataFrame(result,columns=['expression'])
    return result


@app2.callback(Output('umap','figure'),Input('cluster_name','value'))
def update_figure1(cluster_name):
  tempDf = pd.DataFrame()
  tempDf2 = pd.DataFrame()
  if cluster_name == []:
    tempDf = pd.DataFrame(columns=['Cells','UMAP_1','UMAP_2','Cluster'])
    fig = px.scatter(tempDf, x="UMAP_1", y="UMAP_2")
  elif 'All' not in cluster_name:
    for i in cluster_name:
        tempdata1 = umap_data[(umap_data['Cluster'] == i)]
        tempdata2 = umap_data[(umap_data['Cluster'] != i)]
        #print(tempdata1)
        tempDf = tempDf.append(pd.DataFrame(data=tempdata1))
        tempDf1 = tempDf2.append(pd.DataFrame(data=tempdata2))
        tempDf1 = tempDf1.drop(['Cluster'],axis=1)
        tempDf1['Cluster'] = 'other cluster'
        tempDf = pd.concat([tempDf,tempDf1])
        tempDf = tempDf.drop_duplicates()
        #print(tempDf)
        fig = px.scatter(tempDf, x="UMAP_1", y="UMAP_2", color = 'Cluster',color_discrete_map = {'other cluster':'silver'})
        fig.update_traces(marker=dict(size=2))
        fig.update_layout(legend={'itemsizing':'constant'})
  else: 
    tempDf = umap_data
    fig = px.scatter(tempDf, x="UMAP_1", y="UMAP_2", color="Cluster")
    fig.update_traces(marker=dict(size=2))
    fig.update_layout(legend={'itemsizing':'constant'})
  #print(cluster_name)



  #fig = px.scatter(tempDf, x="UMAP_1", y="UMAP_2", color="Cluster")
  return fig


@app2.callback(Output('gene_featureplot','figure'),Input('gene_name','value'))
def update_figure2(gene_name):
  expression = get_expression(gene_name)
  expression_data = pd.concat([umap_data, expression], axis=1)
  expression_data['expression'] = expression_data['expression'].astype(float)
  fig = px.scatter(expression_data, x="UMAP_1", y="UMAP_2", color="expression")
  fig.update_traces(marker=dict(size=2))
  return fig

@app2.callback(Output('TE_featureplot','figure'),Input('TE_name','value'))
def update_figure3(TE_name):
  expression = get_teexpression(TE_name)
  expression_data = pd.concat([umap_data, expression], axis=1)
  expression_data['expression'] = expression_data['expression'].astype(float)
  fig = px.scatter(expression_data, x="UMAP_1", y="UMAP_2", color="expression")
  fig.update_traces(marker=dict(size=2))
  return fig


@app2.callback(Output('linechart1','figure'),Input('gene_name','value'))
def update_figure4(gene_name):
  expression = get_expression(gene_name)
  expression_data = pd.concat([umap_data, expression], axis=1)
  expression_data['expression'] = expression_data['expression'].astype(float)
  y=[]
  for j in range(len(cluster3)):
    for i in cluster3:
      data_j = expression_data[(expression_data['Cluster'] == i)]
      print(data_j)
      #print(data_1)
      mean_j = np.mean(list(data_j['expression']))
      y=y+[mean_j]
    #print(y)
    trace = go.Scatter(x= cluster3, y= y)
    fig = dict(
       data=[trace],
       layout = go.Layout(title= 'changes of gene expression level with stage',
                        yaxis = dict (tickfont = dict(size=13)),
                        xaxis = dict(tickfont = dict (size=10)),
                        margin = dict(l = 110, b = 170, t = 60, r = 30)))
    return fig



@app2.callback(Output('linechart2','figure'),Input('TE_name','value'))
def update_figure5(TE_name):
  expression = get_teexpression(TE_name)
  expression_data = pd.concat([umap_data, expression], axis=1)
  expression_data['expression'] = expression_data['expression'].astype(float)
  y=[]
  for j in range(len(cluster3)):
    for i in cluster3:
      data_j = expression_data[(expression_data['Cluster'] == i)]
      #print(data_j)
      #print(data_1)
      mean_j = np.mean(list(data_j['expression']))
      y=y+[mean_j]
    #print(y)
    trace = go.Scatter(x= cluster3, y= y)
    fig = dict(
       data=[trace],
       layout = go.Layout(title= 'changes of TE expression level with stage',
                        yaxis = dict (tickfont = dict(size=13)),
                        xaxis = dict(tickfont = dict (size=10)),
                        margin = dict(l = 110, b = 170, t = 60, r = 30)))
    return fig


app1lication = DispatcherMiddleware(flaskApp, {'/hPGCLCdb': app1})

if __name__ == '__main__':
    run_simple(hostname="127.0.0.1", port=5050, application=flaskApp)
    # run_server(host="127.0.0.1", port=5050, debug=True)
    # flaskApp.run(host="127.0.0.1", port=5050)
