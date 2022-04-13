import numpy as np
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import math
from pathlib import Path
from scipy.optimize import brentq
import streamlit as st

def read_markdown_file(markdown_file):
    return Path(markdown_file).read_text(encoding="utf8")


@st.cache
def convert_df(df):
   return df.to_csv().encode('utf-8')


class hydro_lsh:
    def __init__(self,Q,loi_frottement='Régime critique',binf_calc=5.,bsup_calc=90.,step=1.0,Lshinf=15.0,Lshsup=40.0):
        self.Q = Q
        self.loi_frottement = loi_frottement
        self.step = step
        self.Lshinf = Lshinf
        self.Lshsup = Lshsup
        self.b_range = np.arange(binf_calc,bsup_calc,step)


        # 'Régime critique','VPE Ferguson (2007)'
        # self.binf = self.Lshinf*self.h_critique_mu(self.Lshinf)
        # self.bsup = self.Lshsup*self.h_critique_mu(self.Lshsup)


    @property
    def binf(self):
        if self.loi_frottement=='Régime critique':
            binf = self.Lshinf*self.h_critique_mu(self.Q,self.Lshinf)
        if self.loi_frottement=='VPE Ferguson (2007)':
            binf = self.Lshinf*self.h_ferguson_mu(self.Q,self.S,self.Lshinf,self.d84)
        return binf

    @property
    def bsup(self):
        if self.loi_frottement=='Régime critique':
            bsup = self.Lshsup*self.h_critique_mu(self.Q,self.Lshsup)
        if self.loi_frottement=='VPE Ferguson (2007)':
            bsup = self.Lshsup*self.h_ferguson_mu(self.Q,self.S,self.Lshsup,self.d84)
        return bsup

    @property
    def h(self):
        if self.loi_frottement=='Régime critique':
            h = self.h_critique_mu(self.Q,self.bsh)
        else:
            h =self.h_ferguson_mu(self.Q,self.S,self.bsh,self.d84)
        return h

    @property
    def b(self):
        return self.bsh*self.h

    @property
    def u(self):
        return self.Q/self.b/self.h
    
    @property
    def H(self):
        return self.h + self.u**2/2/9.81
    
    @property
    def h_range(self):
        if self.loi_frottement=='Régime critique':
            h_range=(self.Q/self.b_range/9.81**0.5)**(2./3.)
        else:
            h_range=np.array([self.h_ferguson(self.Q,self.S,i,self.d84) for i in self.b_range])
        return h_range

    @property
    def bsh_range(self):
        return self.b_range/self.h_range

    @property
    def u_range(self):
        return np.array([self.Q/i/j for i,j in zip(self.h_range,self.b_range)])

    @property
    def H_range(self):
        return self.h_range+self.u_range**2/2/9.81

    @staticmethod
    def h_critique_mu(Q,mu):
        return (Q/mu/9.81**0.5)**(2/5)
    
    @staticmethod
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
    
    @staticmethod
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

    @property
    def Hs_critique(self,b=None):
        if b !=None:
            self.b = b
            self.h = (self.Q/self.b/9.81**0.5)**(2./3.)
            self.A = self.b*self.h
            self.u = self.Q / self.A
            return self.h+(self.Q/self.A)**2/2/9.81
    
    

    
    def make_fig(self,data_trace,rappel=False):
        """make_fig Function that plot Figure

        Args:
            data_trace (list): List containing data series to plot
            rappel (bool, optional): If True plot lines on the plot for specified L/H ratios. Defaults to False.

        Returns:
            figure: Figure
        """

        dict_label={
        'hauteur':{'label':'Hauteur','unit' : '(m)','value':max(self.h_range),'serie':self.h_range},
        'vitesse':{'label':'Vitesse','unit' : '(m/s)','value':max(self.u_range),'serie':self.u_range},
        'charge':{'label':'Charge','unit' : '(m)','value':max(self.H_range),'serie':self.H_range},
        }


        matplotlib.rcParams.update({'font.size':16,'figure.autolayout': True})
        fig=plt.figure(figsize=(12,8))
        gs=gridspec.GridSpec(1,1)
        ax=fig.add_subplot(gs[0])
        ax.grid('--',alpha=0.3)
        
        ### Graph series principales
        if any(i in ['hauteur','vitesse','charge'] for i in data_trace):
            if 'hauteur' in data_trace:
                l1=ax.plot(self.b_range,self.h_range,label="Hauteur")
            if 'vitesse' in data_trace:
                l2=ax.plot(self.b_range,self.u_range,linestyle='dotted',label="Vitesse")
            if 'charge' in data_trace:
                l3=ax.plot(self.b_range,self.H_range,label="Charge")
                
        ### plot labels des axes
            ax.set_xlabel('Largeur du lit (m)')
            data_ax_x = [i for i in data_trace if i !='L sur h']
            s=[dict_label[i]['label'] +' '+dict_label[i]['unit'] for i in data_ax_x]
            s=" - ".join(s)
            ax.set_ylabel(s)
            
        ### Limite axe y
            maxy = max(dict_label[i]['value'] for i in data_ax_x)+1
            ax.set_ylim(0,math.ceil(maxy))

        ### plot lignes de rappel
            if rappel==True:
                ax.plot([self.binf,self.binf],[0,50],'-.r',linewidth=2.0)
                ax.plot([self.bsup,self.bsup],[0,50],'-.r',linewidth=2.0)
                dec_y=0.5
                posy_lsh = (
                    max(
                        np.interp(self.Lshinf, self.bsh_range, j)
                        for j in [dict_label[i]['serie'] for i in data_ax_x]
                    )
                    + dec_y
                )
            
                ax.text(self.binf, posy_lsh,f'L/h={self.Lshinf}' , horizontalalignment='right',rotation=90,fontsize=12)
                ax.text(self.bsup, posy_lsh,f'L/h={self.Lshsup}' , horizontalalignment='right',rotation=90,fontsize=12)
                ax.text(self.binf, 0.1,f'L={round(self.binf,1)} m' , horizontalalignment='right',rotation=90,fontsize=12)
                ax.text(self.bsup, 0.1,f'L={round(self.bsup,1)} m' , horizontalalignment='right',rotation=90,fontsize=12)
                maxy = max(max(dict_label[i]['value'] for i in data_ax_x),posy_lsh)+1
                ax.set_ylim(0,math.ceil(maxy))

        ### recupération des lignes axe principal
        l = [ax.get_lines()[i] for i in range(len(ax.get_lines()))]
        
        ### plot Lsh axe secondaire
        if 'L sur h' in data_trace:
            ax1=ax.twinx()
            l4=ax1.plot(self.b_range,self.bsh_range,'--',color='indigo',label='Rapport L/h',zorder=-1)
            ax1.set_ylabel('Rapport Largeur/hauteur écoulement (-)')
            l = [ax.get_lines()[i] for i in range(len(ax.get_lines()))]+[ax1.get_lines()[i] for i in range(len(ax1.get_lines()))]

        ### plot légende
        labs = [k.get_label() for k in l]
        l=[i for i,j in zip(l,labs) if not j.startswith('_child')]
        labs=[i for i in labs if not i.startswith('_child')]
        leg=ax.legend(l, labs, facecolor='white', framealpha=1,loc=2,fontsize=14,ncol=2)
        
        ### ANNOTATIONS
        annotation_string = r'Q=%.f m$^{3}$/s' % (self.Q)
        ax.text(0.85,0.925,annotation_string,
                horizontalalignment='right',style='normal',
                transform=ax.transAxes,zorder=3,
                bbox={'facecolor':'white', 'alpha':1., 'pad':10},fontsize=14)

        return fig

    def set_d84(self,d84):
        self.d84 = d84

    def set_bsh(self,mu):
        self.bsh = mu  

    def set_S(self,S):
        self.S = S   


        



        