# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 14:44:20 2023

@author: ZimmermannP

"""
import math
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

from src.Refrigerant import refrigerant 
from src.model_steiner import eps
R = 'R717'
ref = refrigerant(R)
g = 9.81


""" Druckverlust nach Chrisholm VDI """

# Reibungsdruckverlust
def deltap_r_x(x,te,m,di):    
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    Re_l0 = m*di/eta_l
    ''' Reibungsdruckverlust nach Chrisholm (VDI) L2 an der Stelle x'''
    # print(self.Re_l0)
    # print(self.Re_g0)
    if Re_l0 < 1187:
        xi_l0 = 64/Re_l0
    else:
        xi_l0 = 0.86859* math.log(Re_l0/(1.964*math.log(Re_l0)-3.825))**-2
    dpdl_l0 = xi_l0 * (m**2/(2*di*rho_l))    
    G = (rho_l/rho_g)*(eta_g/eta_l)**0.2
    phi_l0 = 1 + (G-1) *  ((21/math.sqrt(G))*x**0.9*(1-x)**0.9+x**1.8)
    dpdl = dpdl_l0 * phi_l0
    #print(dpdl)
    return dpdl

def deltap_r_m(m,te,di,t,):
    R = 'R717'
    ref = refrigerant(R)    
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    x1 = 10**-5
    dx = (1-x1)/t
    val = np.empty(shape=(t,2))
    for i in range(t):
        x = x1+dx*i
        dp_r_x = deltap_r_x(x, te, m, di) 
        
        val[i,0] = x
        val[i,1] = dp_r_x

    return val

#Beschleunigungsdruckverlust       
def deltap_b(m,te):
    rho_g = PropsSI('D', 'Q',1,'T',273+te,R) #Dichte gas  [kg/m³]    
    return m**2/rho_g
    
           
def deltap_s_x(theta,x,m):
    ''' geodätische Druckabfall / Länge an der Stelle x '''
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    if theta == 0:
        return 0
    else:
        epsi = eps(x,rho_g, rho_l, m, sigma_r)
        return ((rho_l*(1-epsi) + rho_g*epsi)*g*math.sin(math.radians(theta))) 
       
def deltap_r_G(te,m,di,x1,x2=1,):
    """ Druckabfall nach Gronerund Integralmittelwert"""
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    dpdl = 0.158 * ((eta_l**0.25*m**1.75)/(di**1.25*rho_l)) #flüssiger druckverlust
    
    #2Phasen-Multiplikator
    O = (rho_l/rho_g) / (eta_l/eta_g)**0.25 -1
    Fr_l = m**2 / (rho_l**2*di*g)
    if Fr_l >= 1: ffr = 1
    else: ffr = Fr_l**0.3 + 0.0055 * (math.log(1/Fr_l))**2
    #Integralmittelwert
    dpdl_fr = (ffr/(1-x1)) * ((0.5*(1-x1**2)) + 1.429 * (1-x1**2.8) - 0.364* ffr**0.5 * (1-x1**11))
    Ph = O*dpdl_fr+1
    
    dpdl_tf = dpdl * Ph
    return dpdl_tf

def deltap_r_L(te,m,di,x1,L,p,x2=1):
    """ mittlerer Druckabfall * Rohrlänge * Pässe """
    dpL = deltap_r_G(te,m,di,x1,x2=1,)
    dp = dpL * L * p
    return dp
    
    
    
    
    
if __name__ == "__main__":
    R = 'R717'
    
    x = 0.9
    te = -30
    m = 15 
    di = 0.015
    t = 100
    dp_r = deltap_r_x(x, te, m, di)
    
    val = deltap_r_m(m, te, di, t)
    xd = val[:,0]
    dpr1 = val[:,1] / 100
    
    plt.plot(xd,dpr1, label = (f'Stahl {di*1000:.2f} x {0.0005*1000}mm Glattrohr'))
    plt.xlim(0,1)
    plt.ylim(0,max(dpr1)*1.1)
    plt.ylabel('\u0394p / \u0394l [mbar/m]')
    plt.xlabel('Massendampfgehalt [-]')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
    plt.title(f'{R};  te={te}°C; Glattrohr')
    plt.grid(True)
