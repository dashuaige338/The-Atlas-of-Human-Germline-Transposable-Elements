import os
os.chdir('F:\\LabW')
import dash
app = dash.Dash(__name__, suppress_callback_exceptions=True)
app.title = 'NANOS3 Specific Gene'
server = app.server