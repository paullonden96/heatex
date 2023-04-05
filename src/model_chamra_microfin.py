# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:56:12 2023

@author: ZimmermannP
"""

""" Modell für microfin-tubes nach Chamra 
    Quelle: https://journals.sagepub.com/doi/epdf/10.1243/0954406JMES131 """
from CoolProp.CoolProp import PropsSI
import math
import matplotlib.pyplot as plt
#from mpmath import coth
import numpy as np
#import pandas as pd
import time
#from sympy import solve,solveset, Symbol, S, sin, Eq
import random
from src.Refrigerant import refrigerant 


Pi = math.pi
g = 9.81

class abschnitt():
    
    def __init__(self,R:str,te,x,m,di,s,theta,lambda_w,q,lambda_g,lambda_l,
                 rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,y,e,nf,b,th):
        
        self.R = R # Kältemittel
        self.te = te # Verdampfungstemperatur
        self.x = x # Dampfgehalt im Abschnitt
        self.m = m # Massenstromdichte
        self.di = di # Rohrinnendurchmesser
        self.theta = theta # Neigungswinkel zur Horizontalen
        self.q = q # innere Wärmestromdichte
        self.s = s # Wandstärke
        self.lambda_w = lambda_w # Wärmeleitfähigkeit 
        self.lambda_g = lambda_g # Wärmeleitfähigkeit Gas
        self.lambda_l = lambda_l # Wärmeleeitfähigkeit Flüssigkeit
        self.rho_l = rho_l # Dichte Flüssigkeit
        self.rho_g = rho_g # Dichte Gas
        self.eta_g = eta_g # Viskosität Gas
        self.eta_l = eta_l # Viskosität Flüssigkeit
        self.pr_g = pr_g # Prandtl Gas
        self.sigma_r = sigma_r # Oberflächenspannung
        self.cp_g = cp_g # spez. Wärmekapazität Gas
        self.cp_l = cp_l # spez. Wärmekapazität Flüssigkeit
        self.delta_hv = delta_hv # Enthalpiedifferenz Verdampfen 
        self.pe = pe # Verdampfungsdruck
       
        # Parameter von microfin tubes 
        self.y = y # helix angle gamma
        self.e = e # fin height
        self.nf = nf # number of microfins per unit lenght
        self.b = b # angle beta apex angle
        self.th = th # tube thickness
        
       
        
        self.Re_l0 = self.m*self.di/self.eta_l
        self.Re_g0 = self.m*self.di/self.eta_g
        self.Re_l = self.m*self.x*self.di/eta_g
        self.Re_g = self.m*(1-self.x)*self.di/self.eta_l
        
        self.Pr_l = PropsSI('Prandtl','Q',0,'T',te+273.15,R)
        self.Pr_g = PropsSI('Prandtl','Q',1,'T',te+273.15,R)
        
        # Martinelli Parameter --> eventuell 0.875 --> 0.9
        # 0.125 --> 0.1 
        self.X = ((1-self.x)/self.x)**0.9 * (self.rho_g/self.rho_l)**0.5 *\
            (self.eta_l / self.eta_g)**0.1
                
    def alpha_tf_x(self):
        #nucelar boiling --> prüfen ob überhaupt vorliegt?
        q = self.delta_hv*self.m
        h_pb = h_nb(self.pe,self.q)
        
        #h_pb = 0
        # flow boiling effect factor
        F1 = 0.01*self.di**-1
        
        # liquidphase heat transfer ?? why Dittus Boiter --> vielleicht von zürcher nehmen
        C = 0.01361
        m = 0.6965
        Re_s =  (self.m*(1-self.x)*self.di) / self.eta_l
        hl = 4*C*Re_s**m*self.Pr_l**0.4*(self.lambda_l/self.di)
        print(hl)
        
        hl2 = 0.023*Re_s**0.8*self.Pr_g**0.4*(self.lambda_l/self.di)
        print(hl2)
        #two-phase multiplier
        Phi = ((1-self.x)+2.63*self.x*(self.rho_l/self.rho_g)**0.5)**0.8
        #Geometry factor
        Rx = 1/math.cos(self.y) * (((2*self.e*self.nf*(((1-math.sin(self.b/2)))))/(Pi*self.di*math.cos(self.b/2)))+1)
        #Bond Number
        Bon=(g*self.rho_l*self.e*Pi*self.di)/(8*self.sigma_r*self.nf)
        #Froud Number
        Frv = self.m**2 /(self.rho_g**2*g*self.di)
        #flow parameters 
        F2 = 0.01*self.di**-1
        F3 = 100*self.m**-1
        
        C1 = 1.516
        C2 = 1.161
        C3 = -1.7640
        C4 = 2.662
        C5 = -0.2158
        C6 = 0.5927  
        C7 = 0.0582
        
        print('h_pb=',h_pb)
        print('Rx =',Rx)
        return h_pb * (C1*self.X**C2)*F1**C3 +\
            hl2 * Phi * Rx**C4 * (Bon * Frv)**C5 * F2**C6 * F3**C7
                
def h_nb(pe,ql):
    """ nuclear boiling Cooper """
    p_c = 113.33 * 10**5
    M0 = 17.03
    return 55 * (pe/p_c)**0.12 * (-1*math.log(pe/p_c))**-0.55*M0**-0.5*ql**0.67

def Kabelac():
    """ Parameter aus den Experimenten von microfin mit Ammoniak von Kabelac """
    # square shaped --> beta = 0?
    nf = 21 # Anzahl an Rippen pro Umfang
    y = 25 # Drallwinkel Gamma
    e = 0.63/1000 # Zahnhöhe 
    b = 0 # Winkel einr Rippe
    th = 0.1/1000 # Wanddicke Rohr von der Zahnwurzel gemessen
    di = 11.13 / 1000 # Innendurchmesser
    return y,e,nf,b,th

def Zürcher():
    """ Parameter aus den Experimenten von microfin mit Ammoniak von Zürcher """
    nf = 34
    y = 18
    e = 0.33/1000
    b = 0
    th = 0.1/1000
    #di = 0.01346
    return y,e,nf,b,th

def Kelly():
    # di = 0.0118
    y = 17
    nf = 60
    e = 0.25/1000
    th = None
    b = None
    return y,e,nf,b,th
    
    

def alpah_tf_m(m,q,di,s,te,t,geo,R='R717', theta=0, lw=26):
    ref = refrigerant(R)
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R) 
    pe  = PropsSI('P','T',273.15+te,'Q',1,R) #Verdampfungsdruck
    x1 = 0.1
    dx = (1-x1)/t
    val = np.empty(shape=(t,2))
    for i in range(t):
        x = x1+(dx*i)
        print(f'Massendampfgehalt = {x}')
        a = abschnitt(R,te,x,m,di,s,theta,lw,q,lambda_g,lambda_l,
                      rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,*Kabelac())
        alpha_tf = a.alpha_tf_x()
        val[i,0] = x
        val[i,1] = alpha_tf
        
    return val
          
if __name__ == "__main__":
    R = 'R717'
    m = 50
    di = 0.01113
    s = 0.5 / 1000
    te = -20
    Ra = 35 * 10**-6
    lw = 16 # Wärmeleitfähigkeit Rohr
    
    
    x = 0.5
    theta = 0 #Winkel zur Horizontalen
    pe  = PropsSI('P','T',273+te,'Q',1,R) #Verdampfungsdruck
    # q = m*delta_hv*delta_x?
    t = 20
    
    alpha_tf_x()
    geo = Kabelac()
    val = alpah_tf_m(m, di, te, t, geo)
    
    xd = val[:,0]
    alpha_z = val[:,1]
    xd2 = [0.16,0.34,0.52,0.9]
    alpha_z2 = [12000,14500,16000,13000]
    plt.plot(val[:,0], val[:,1], label = (f'berechnet'))
    plt.scatter(xd2,alpha_z2,label = (f'Experimentielle Werte nach Kabelac'))
    plt.ylim(0,max(alpha_z2)*1.1)
    plt.xlim(0,1)
    plt.ylabel('\u03B1_i [W/m²K]')
    plt.xlabel('Massendampfgehalt [-]')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
    #plt.title(f'{R}; {te}°C;')
    plt.title(f'{R}; te= {te}°C; {m} [kg/m²s] Rippenrohr, Modell nach Chamra, Rippenrohr,')
    plt.grid(True)
    
    
            
            
            