import trialObject
import numpy as np
import plotly.graph_objects as go

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input,Output,State

external_stylesheets = [{'href':'https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css',
                        'rel':'stylesheet'}]
app = dash.Dash(__name__,title='Value-based clinical trials',external_stylesheets=external_stylesheets)
server=app.server



trial = trialObject.Trial()

with open('manual.md','r') as f:
    manual = f.read()

with open('intro.md','r') as f:
    intro = f.read()

app.layout = html.Div(className='app-layout', children=[
    html.Div(className='container',children=[
        html.H1('Value-based clinical trials'),
        dcc.Markdown(intro),
        html.Div(className = 'row', children=[
            html.Div(className = 'col-md-6 bg1', children=[
                html.H3('Parameters'),
                html.Div(className = 'col-md-6', children=[
                    html.H4('Population'),
                    html.Div(id='population-formula'),
                    html.Div(className='input',children=[dcc.RadioItems(id='FixedPool',options=[{'label':'Fixed pool P=ζH','value':1},{'label':'Fixed horizon  P(T)=ζ(H-T-Δ)','value':0}],value=int(trial.FixedPool),labelStyle={'paddingRight':'10px', 'fontWeight':'normal'})]),
                    html.Div(className='input',children=[html.Div('Time horizon H (months)'),dcc.Input(id='H',type='number',min=0,value=trial.H)]),
                    html.Div(className='input',children=[html.Div('Incidence rate ζ (per month)'),dcc.Input(id='incidence',type='number',min=0,value=trial.incidence)]),
                    # html.Div(className='input',children=['Incidence rate ζ = ',html.Span(id='incidencevalue'),dcc.Slider(id='incidence',min=10.0,max=1000.0,value=trial.incidence,step=10,marks={i:str(i) for i in [10,1000,200,400,600,800]})]),
                    html.Div(className='input',children=[html.Div(['Delay Δ (months)']),dcc.Input(id='delay',type='number',value=trial.delay)]),
                    # html.Div(className='input',children=['Discount rate ρ = ',html.Span(id='discratevalue'),'%',dcc.Slider(id='discrate',min=0,max=0.3,value=trial.discrate,step=0.005,marks={i:str(i*100) for i in [0,0.1,0.2,0.3]})]),
                    html.Div(className='input',children=['p',html.Sub('N'),'=',html.Span(id='pNvalue'), dcc.Slider(id='pN',min=0,max=1,value=trial.pN,step=0.01,marks={i:'{:.1f}'.format(i) for i in [0,0.2,0.4,0.6,0.8,1]})]),
                    html.Div(className='input',children=[dcc.RadioItems(id='online',options=[{'label':'offline learning','value':0},{'label':'online learning','value':1}],value=trial.online,labelStyle={'paddingRight':'10px', 'fontWeight':'normal'})]),
                    html.Div(className='input',children=[html.Div('Discount rate ρ (per month)'),dcc.Input(id='discrate',type='number',min=0,value=trial.discrate)]),
                ]),
                html.Div(className = 'col-md-6 left-line', children=[
                    html.H4('Costs'),
                    html.Div(className='input',children=[html.Div('c (sampling cost)'),dcc.Input(id='c',type='number',min=0,value=trial.c)]),
                    html.Div(className='input',children=[html.Div(['Cost of capacity c',html.Br(),html.Sub('cap'),'(r)= c',html.Sub('r'),'*r',html.Sup('θ'),'=',html.Span(id='crate-value'),'*r',html.Sup(id='theta-value')])]),
                    html.Div(className='input',children=[html.Div(['c',html.Sub('r'),'=']),dcc.Input(id='crate',type='number',min=0,value=100)]),
                    html.Div(className='input',children=[html.Div('θ='),dcc.Input(id='theta',type='number',min=1,value=2)]),
                    html.Div(className='input',children=[html.Div(['I',html.Sub('N'), ' (Cost of switching to N)']),dcc.Input(id='IN',type='number',value=trial.IN)]),
                    html.Div(className='input',children=[html.Div(['I',html.Sub('S'), ' (Cost of switching to S)']),dcc.Input(id='IS',type='number',value=trial.IS)]),
                    html.H4('Statistical parameters'),
                    html.Div(className='input',children=[html.Div(['μ',html.Sub('0'), ' (prior mean)']),dcc.Input(id='mu0',type='number',value=trial.mu0)]),
                    html.Div(className='input',children=[html.Div(['σ',html.Sub('X'), ' (sample standard deviation)']),dcc.Input(id='sigma',type='number',min=0,value=trial.sigma)]),
                    html.Div(className='input',children=[html.Div(['n',html.Sub('0'), ' (prior sample size)']),dcc.Input(id='n0',type='number',value=trial.n0)]),
                ]),
            ]),
            html.Div(className = 'col-md-6', children=[
                dcc.Tabs(className='custom-tabs',value='tab-1',children=[
                    dcc.Tab(className='custom-tab',label='Plot',value='tab-1',children=[
                        html.H3('Expected net gain V(T,r)'),
                        dcc.Graph(id='VOI-figure'),
                        html.Div(className='input',children=[dcc.RadioItems(id='surf',options=[{'label':'Contour plot','value':0},{'label':'Surface plot','value':1,'disabled':True}],value=0,labelStyle={'paddingRight':'10px', 'fontWeight':'normal'})]),
                        html.Button(id='maximizer',className='btn btn-primary btn-block',children='Show optimal design'),
                        html.Div(className = 'row',children=[
                            html.Div(className = 'col-md-3', children=['V*=',html.Span(id='Vopt')]),
                            html.Div(className = 'col-md-3', children=['T*=',html.Span(id='Topt')]),
                            html.Div(className = 'col-md-3', children=['r*=',html.Span(id='ropt')]),
                            html.Div(className = 'col-md-3', children=['Q*=',html.Span(id='Qopt')]),
                        ]),
                        html.Div(className = 'row',children=[
                            html.Div(className = 'col-md-6', children=[html.Div(className='input',children=[html.Div(['T',html.Sub('max')]),dcc.Input(id='Tmax',type='number',min=0,value=trial.Tmax)])]),
                            html.Div(className = 'col-md-6', children=[html.Div(className='input',children=[html.Div(['r',html.Sub('max')]),dcc.Input(id='rmax',type='number',min=0,value=trial.rmax)])])
                        ]),
                    ]),
                    dcc.Tab(className='custom-tab',label='Manual',value='tab-2',children=[
                        dcc.Markdown(manual)
                    ])
                ]),
            ]),
        ])
    ]),
    html.Div(style={'height': '40px'}) # White space at the bottom
])



@app.callback([Output('VOI-figure', 'figure'),Output('Vopt','children'),Output('Topt','children'),
            Output('ropt','children'),Output('Qopt','children')],
            [Input('maximizer', 'n_clicks'),
            Input('incidence', 'value'), Input('c','value'),Input('pN','value'),Input('H','value'),
            Input('Tmax','value'),Input('rmax','value'),Input('discrate','value'),
            Input('sigma','value'),Input('mu0','value'),Input('n0','value'),Input('online','value'),
            Input('crate','value'),Input('theta','value'),Input('delay','value'),Input('IN','value'),Input('IS','value'),
            Input('FixedPool','value'),
            Input('surf','value')])
def updateFigure(n_clicks,incidence,c,pN,H,Tmax,rmax,discrate,sigma,mu0,n0,online,crate,theta,delay,IN,IS,FixedPool,surf):
    trial.c = c
    trial.H = H
    trial.incidence = incidence
    trial.pN = pN
    trial.Tmax = Tmax
    trial.rmax = rmax
    trial.discrate = discrate
    trial.sigma = sigma
    trial.mu0 = mu0
    trial.n0 = n0
    trial.online = online
    trial.ccap = lambda r: crate*r**theta
    trial.delay = delay
    trial.IN = IN
    trial.IS = IS
    trial.FixedPool = FixedPool

    user_click = dash.callback_context.triggered[0]['prop_id'].split('.')[0]
    if not user_click or user_click != 'maximizer':
        return trial.plot(max=False,surf=surf),'','','',''
    else:
        fig = trial.plot(max=False,surf=surf)
        res=trial.maximize()
        fig.add_trace(go.Scatter(x=[res.x[0]], y=[res.x[1]], name=str(trial.value(res.x[0],res.x[1])),marker_color='black',marker_size=10))
        return fig,'{:.3f}M'.format(trial.value(res.x[0],res.x[1])/1e6),'{:.2f}'.format(res.x[0]),'{:.2f}'.format(res.x[1]),'{:.2f}'.format(res.x[0]*res.x[1]/2)

@app.callback(Output('population-formula','children'),
            Input('FixedPool','value'),Input('H','value'),
            Input('incidence','value'),Input('delay','value'))
def updatePopFormula(FixedPool,H,incidence,delay):
    if FixedPool:
        return 'P=' + str(incidence*H)
    else:
        return 'P(T)={}*({}-T)'.format(incidence,H-delay)

@app.callback(Output('theta-value','children'),[Input('theta','value')])
def updatetheta(x):
    return str(x)

@app.callback(Output('crate-value','children'),[Input('crate','value')])
def updatecrate(x):
    return str(x)

# @app.callback(Output('incidencevalue','children'),[Input('incidence','value')])
# def updateIncidence(x):
#     return str(x)

@app.callback(Output('pNvalue','children'),[Input('pN','value')])
def updatepN(x):
    return str(x)

# @app.callback(Output('discratevalue','children'),[Input('discrate','value')])
# def updateDiscrate(x):
#     return str(x*100)



if __name__ == '__main__':
    app.run_server(debug=True)