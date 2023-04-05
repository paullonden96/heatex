# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 15:38:50 2023

@author: ZimmermannP
"""

""" microfin Korrelaion nach Thome
    https://www.researchgate.net/publication/37421421_Evaporation_in_Microfin_Tubes_A_Generalized_Prediction_Model?enrichId=rgreq-376250fb3a800b22a0cb9a73e48e789c-XXX&enrichSource=Y292ZXJQYWdlOzM3NDIxNDIxO0FTOjIyODAyNjUyNDg5MzE4NUAxNDMxMzc3MTY4NDgx&el=1_x_2&_esc=publicationCoverPdf """

from CoolProp.CoolProp import PropsSI
import math
import matplotlib.pyplot as plt
#from mpmath import coth
import numpy as np
#import pandas as pd
#import time
#from sympy import solve,solveset, Symbol, S, sin, Eq
#import random
from src.Refrigerant import refrigerant 
from src.model_steiner import eps
from src.model_chamra_microfin import Kabelac

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
        self.Re_g = self.m*self.x*self.di/eta_g
        self.Re_l = self.m*(1-self.x)*self.di/self.eta_l
        
        self.Pr_l = PropsSI('Prandtl','Q',0,'T',te+273.15,R)
        self.Pr_g = PropsSI('Prandtl','Q',1,'T',te+273.15,R)
        
        # Martinelli Parameter --> eventuell 0.875 --> 0.9
        # 0.125 --> 0.1 
        self.X = ((1-self.x)/self.x)**0.9 * (self.rho_g/self.rho_l)**0.5 *\
            (self.eta_l / self.eta_g)**0.1
                
    def alpha_tf_x(self):
        C = 0.01361
        m = 0.6965
        #enhancement factor ripped tubes
        pf =  (Pi*self.di) / self.nf #axial pitch from fin to fin 
        
        E_rb = (1+(2.64*self.Re_l**0.036*(self.e/self.di)**0.212*(pf/self.di)**-0.21*\
                  (self.y/90)**0.29*self.Pr_l**-0.024)**7)**(1/7)
        #print('enhancement=',E_rb)
        epsi = eps(self.x, self.rho_g, self.rho_l, self.m, self.sigma_r)
        sigma = self.di*(1-epsi) / 4 
        Re_s = 4*self.m*(1-self.x)*sigma/((1-epsi)*self.eta_l)
        #h_cb_l = 4*C*Re_s**m*self.Pr_l**0.4*(self.lambda_l/sigma)
        h_cb_l = C*Re_s**m*self.Pr_l**0.4*(self.lambda_l/sigma)
        return E_rb * h_cb_l
    
def alpha_tf_m_thome(m,q,di,s,te,t,geo,R='R717', theta=0, lw=26):
    ref = refrigerant(R)
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R) 
    pe  = PropsSI('P','T',273.15+te,'Q',1,R) #Verdampfungsdruck
    x1 = 0.1
    dx = (1-x1)/t
    val = np.empty(shape=(t,2))
    for i in range(t):
        x = x1+(dx*i)
        #print(f'Massendampfgehalt = {x}')
        a = abschnitt(R,te,x,m,di,s,theta,lw,q,lambda_g,lambda_l,
                      rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,*Kabelac())
        alpha_tf = a.alpha_tf_x()
        val[i,0] = x
        val[i,1] = alpha_tf
        
    return val
if __name__ == "__main__":
    geo = Kabelac()
    R = 'R717'
    m = 50
    te = -20
    q = 40000
    di = 0.01113
    s = 0.5 / 1000
    t = 50
    val = alpha_tf_m_thome(m,q, di,s, te, t, geo)
    
    xd = val[:,0]
    alpha_z = val[:,1]
    xd2 = [0.16,0.34,0.52,0.9]
    alpha_z2 = [12000,14500,16000,13000]
    plt.plot(val[:,0], val[:,1], label = (f'berechnet'))
    plt.scatter(xd2,alpha_z2,label = (f'Experimentielle Werte nach Kabelac'))
    plt.ylim(0,max(alpha_z)*1.1)
    plt.xlim(0,1)
    plt.ylabel('\u03B1_i [W/m²K]')
    plt.xlabel('Massendampfgehalt [-]')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
    #plt.title(f'{R}; {te}°C;')
    plt.title(f'{R}; te= {te}°C; {m} [kg/m²s] Rippenrohr, Modell nach Thome, Rippenrohr,')
    plt.grid(True)