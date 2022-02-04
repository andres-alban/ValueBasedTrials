import numpy as np
from scipy.stats import norm
import plotly.graph_objects as go
from scipy.optimize import minimize,Bounds,BFGS

def psi(x):
    y = norm.pdf(x) - x*(1-norm.cdf(x))
    if y.size==1:
        if np.isnan(y): return 0
    else:
        y[np.isnan(y)] = 0
    return y


class Trial():
    def __init__(self,c=1000,pN=0.25,IN=0,IS=0,H=200,incidence=100,delay=0, ccap=lambda r: 100*r**2,Tmax=100,rmax=100,discrate=0.01/200,mu0=0,n0=1,sigma=1000,online=0,FixedPool=True):
        self.c = c
        self.pN = pN
        self.IN = IN
        self.IS = IS
        self.H = H
        self.incidence = incidence
        self.delay = delay
        self.ccap = ccap
        self.Tmax = Tmax
        self.rmax = rmax
        self.discrate = discrate
        self.mu0 = mu0
        self.n0 = n0
        self.sigma = sigma
        self.online = online
        self.FixedPool = FixedPool


    def sigmaZ(self,Q):
        return np.sqrt(self.sigma**2 * Q / (self.n0 * (self.n0 + Q)))

    def Teff(self,T):
        if self.discrate == 0:
            return T
        else:
            return (1 - np.exp(-self.discrate*T))/self.discrate

    def Pop(self,T):
        if self.FixedPool:
            return self.incidence*self.H
        else:
            return self.incidence*(self.H - self.delay - T)

    def value(self,T,r,mesh=False):
        if mesh:
            TT, rr = np.meshgrid(T,r)
        else:
            TT=T
            rr=r

        if self.discrate>0:
            P = self.incidence/self.discrate*(1-np.exp(-self.discrate*self.Pop(TT)/self.incidence))
        else:
            P=self.Pop(TT)

        if self.IN == 0:
            alphaN = 0
        elif self.pN == 0:
            alphaN = np.inf
        else:
            alphaN = self.IN/((1-self.pN)*P)

        if self.IS == 0:
            alphaS = 0
        elif self.pN == 1:
            alphaS = np.inf
        else:
            alphaS = self.IS/(self.pN*P)


        y = ( 0.5*self.online*rr*self.Teff(TT)*(1-2*self.pN)*self.mu0
        + np.exp(-self.discrate*(TT + self.delay))*P*(1 - self.pN)*self.sigmaZ(TT*rr/2)*psi((alphaN - self.mu0)/self.sigmaZ(TT*rr/2))
        + np.exp(-self.discrate*(TT + self.delay))*P*self.pN*self.sigmaZ(TT*rr/2)*psi((alphaS + self.mu0)/self.sigmaZ(TT*rr/2))
        - self.ccap(rr) - self.c*rr*self.Teff(TT) )

        return y

    def plot(self,surf=False,max=True):
        T = np.linspace(0,self.Tmax,100)
        r = np.linspace(0,self.rmax,100)
        z=self.value(T,r,mesh=True)
        fig = go.Figure()
        if surf:
            fig.update_layout(scene=dict(
            zaxis_title='Expected net gain',
            xaxis_title='Trial length (T)',
            yaxis_title='recruitment rate (r)'),
            margin=dict(l=0,r=0),
            paper_bgcolor='white',
            plot_bgcolor='white'
            )
            fig.add_trace(go.Surface(z=z,x=T,y=r,name='Value function'))
            if max == True:
                res=self.maximize()
                fig.add_trace(go.Scatter3d(x=[res.x[0]], y=[res.x[1]], z=[self.value(res.x[0],res.x[1])],name='Optimal value'))
        else:
            fig.update_layout(
            xaxis_title='Trial length (T)',
            yaxis_title='recruitment rate (r)',
            margin={'t':40,'b':10},
            height=300
            )
            fig.add_trace(go.Contour(z=z,x=T,y=r,name='Value function'
            # , contours_coloring='heatmap',contours=dict(start=np.amin(z),end=np.amax(z),size=(np.amax(z)-np.amin(z))/10)
            # , colorscale=[[0, '#fee0d2'],[0.5,'#fc9272'],[1, '#de2d26']]
            , contours=dict(showlabels=True)
            ))
            if max == True:
                res=self.maximize()
                fig.add_trace(go.Scatter(x=[res.x[0]], y=[res.x[1]], name=str(self.value(res.x[0],res.x[1])),marker_color='black',marker_size=10))
        return fig

    def objective(self,x):
        T,r = x
        return -self.value(T,r)

    def maximize(self,T0=1,r0=1):
        bounds = Bounds([0,0],[self.Tmax,self.rmax])
        return minimize(self.objective,[T0,r0],method='trust-constr',jac='2-point',hess=BFGS(),bounds=bounds)

