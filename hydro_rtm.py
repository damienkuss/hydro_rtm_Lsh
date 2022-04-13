from pickle import BINFLOAT
from scipy.optimize import brentq
import numpy as np


def h_ferguson_mu(Q,S,mu,d84):
    def Er(h,args):
        h=h+0.000001
        A=mu*h**2
        P=(mu+2)*h
        R=A/P
        Qstar=mu*h**2*(9.81*S*R)**0.5*2.5*(R/d84)/(1+0.15*(R/d84)**(5/3))**0.5
        return (Qstar-Q)
    if Q==0:
        h=0
    else:
        h=brentq(Er,0.000001,1000.,args=([Q,S,mu,d84]))
    return h

def Hs_ferguson_mu(Q,S,mu,d84):
    def Er(h,args):
        h=h+0.000001
        A=mu*h**2
        P=(mu+2)*h
        R=A/P
        Qstar=mu*h**2*(9.81*S*R)**0.5*2.5*(R/d84)/(1+0.15*(R/d84)**(5/3))**0.5
        return (Qstar-Q)
    if Q==0:
        h=0
        Hs=0
    else:
        h=brentq(Er,0.000001,1000.,args=([Q,S,mu,d84]))
        A=mu*h**2
        Hs=h+(Q*A)**2/2/9.81
    return Hs

def h_ferguson(Q,S,b,d84):
    def Er(h,args):
        h=h+0.000001
        A=b*h
        P=b+2*h
        R=A/P
        Qstar=b*h*(9.81*S*R)**0.5*(2.5*(R/d84))/(1+0.15*(R/d84)**(5/3))**0.5
        return (Qstar-Q)
    if Q==0:
        h=0
    else:
        h=brentq(Er,0.000001,1000.,args=([Q,S,b,d84]))
    return h

def Q_ferguson(h,S,b,d84):
    h=h+0.000001
    A=b*h
    P=b+2*h
    R=A/P
    return b*h*(9.81*S*R)**0.5*2.5*(R/d84)/(1+0.15*(R/d84)**(5/3))**0.5

def h_critique(Q,b):
    return (Q/b/9.81**0.5)**(2./3.)

def Q_critique(h,b):
    return b*9.81**0.5*h**1.5

def Hs_critique(Q,b):
    h=(Q/b/9.81**0.5)**(2./3.)
    A=b*h
    return h+(Q/A)**2/2/9.81

def Hs_ferguson(Q,S,b,d84):
    def Er(h,args):
        h=h+0.000001
        A=b*h
        P=b+2*h
        R=A/P
        Qstar=b*h*(9.81*S*R)**0.5*2.5*(R/d84)/(1+0.15*(R/d84)**(5/3))**0.5
        return (Qstar-Q)
    if Q==0:
        h=0
        Hs=0
    else:
        h=brentq(Er,0.000001,1000.,args=([Q,S,b,d84]))
        A=b*h
        Hs=h+(Q/A)**2/2/9.81
    return Hs

def h_critique_mu(Q,mu):
    return (Q/mu/9.81**0.5)**(2/5)

def Hs_critique_mu(Q,mu):
    h=(Q/mu/9.81**0.5)**(2/5)
    A=mu*h**2
    return h+(Q/A)**2/2/9.81

class hydro_lsh():
    def __init__(self,Q,binf,bsup,step=1.0,Lshinf=15.0,Lshsup=40.0):
        self.binf = binf
        self.bsup = bsup
        self.Q = Q
        self.step = step
        self.Lshinf = Lshinf
        self.Lshsup = Lshsup
        self.b = np.arange(binf,bsup,step)
        self.h=(self.Q/self.b/9.81**0.5)**(2./3.)
        self.bsh=self.b/self.h
        self.u=(9.81*self.h)**0.5
        self.H=self.h+self.u**2/2/9.81
        binf=Lshinf*h_critique_mu(Q,Lshinf)
        bsup=Lshsup*h_critique_mu(Q,Lshsup)

    def h_critique_mu(Q,mu):
        return (Q/mu/9.81**0.5)**(2/5)

    