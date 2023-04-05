# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 14:44:20 2023

@author: ZimmermannP

"""
import math
import numpy as np
import matplotlib.pyplot as plt

from src.Refrigerant import refrigerant 


""" Druckverlust nach Chrisholm VDI """

# Reibungsdruckverlust
def deltap_r_x(x,te,m,di):
    R = 'R717'
    ref = refrigerant(R)    
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    
    Re_l0 = m*di/eta_l
    Re_g0 = m*di/eta_g
    
   
    
    ''' Reibungsdruckverlust nach Chrisholm (VDI) L2 an der Stelle x'''
    # print(self.Re_l0)
    # print(self.Re_g0)
    if Re_l0 or Re_g0 < 1055:
        xi_l0 = 64/Re_l0
        xi_g0 = 64/Re_g0
    else:
        xi_l0 = (0.86859* math.log10(Re_l0/(1.964*math.log10(Re_l0)-3.825)))**-2
        xi_g0 = (0.86859* math.log(Re_g0/(1.964*math.log10(Re_g0)-3.825)))**-2
    dpdl_l0 = xi_l0 * (m**2/(2*di*rho_l))
   # dpdl_g0 = xi_g0 * (m**2/(2*di*rho_g))    
    G = (rho_l/rho_g)*(eta_g/eta_l)**0.2
    phi_l0 = 1 + (G-1) *  ((21/math.sqrt(G))*x**0.9*(1-x)**0.9+x**1.8)
    dpdl = dpdl_l0 * phi_l0
    #print(dpdl)
   
    return dpdl


def deltap_r_m(m,te,di,t):
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
        print(i)
        print(x)
    return val
        
    
    
if __name__ == "__main__":
    x = 0.9
    te = -20
    m = 100 
    di = 0.0081
    t = 100
    dp_r = deltap_r_x(x, te, m, di)
    
    val = deltap_r_m(m, te, di, t)
    xd = val[:,0]
    dpr1 = val[:,1]
    
    plt.plot(xd,dpr1, label = (f'Stahl {di*1000:.2f} x {0.0005*1000}mm Glattrohr'))
    plt.xlim(0,1)
    plt.ylim(0,max(dpr1+200))
    plt.ylabel('\u03B1_i [W/m²K]')
    plt.xlabel('Massendampfgehalt [-]')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
    #plt.title(f'{R}; {te}°C;')
    plt.title(f'{R};  te={te}°C; q={q/1000} kW/m² Glattrohr')
    plt.grid(True)
