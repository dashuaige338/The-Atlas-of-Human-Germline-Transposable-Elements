import os
os.chdir('F:\\LabW')
import dash_core_components as dcc
import dash_html_components as html
import dash
from xzkapp import app
from xzkapp import server
from xzklayout import totals
import xzkcallback
import dash

app.layout = totals


if __name__ == "__main__":
    app.run_server(debug=True)