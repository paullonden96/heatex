# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:09:27 2023

@author: ZimmermannP
"""

""" einfaches ploten von Graphen für thermische Eigenschaften von Ammoniak """

import numpy as np
import matplotlib.pyplot as plt

from src.Refrigerant import refrigerant 
from src.model_steiner import eps



def plot_graph(x,y,l,ylabel,xlabel,title):
    plt.plot(x,y,label = l)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
    plt.grid(True)
    plt.xlim(0)
    plt.ylim(0)
        
    
    
if __name__ == "__main__":

    R = 'R717'
    ref = refrigerant(R)
    
    
    """ Volumetrischer Dampfgehalt in Abhängigkeit der Verdampfungstemperatur nach Rouhani """
    te1 = -40
    te2 = -20
    te3 = 0
    
    t = 100
    m = 15
    dx = (1-0.001)/t
    val = np.empty(shape=(t,4))
    for i in range(t):
        x = 0.001+dx*i
        lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te1, R)
        epsi1 = eps(x, rho_g, rho_l, m, sigma_r)
        lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te2, R)
        epsi2 = eps(x, rho_g, rho_l, m, sigma_r)
        lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te3, R)
        epsi3 = eps(x, rho_g, rho_l, m, sigma_r)
        
        val[i,0] = x
        val[i,1] = epsi1
        val[i,2] = epsi2
        val[i,3] = epsi3
    
    xd = val[:,0]
    plt.plot(xd,val[:,1],label ='-40 °C')
    plt.plot(xd,val[:,2],label='-20 °C')
    plt.plot(xd,val[:,3],label= '0 °C')
    plt.ylabel('Volumetrischer Dampfgehalt \u03B5 [-]')
    plt.xlabel('Massendampfgehalt x [-]')
    plt.title(f'{R}; nach Rouhani (heterogenes Modell); {m} kg/m²s')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.grid(True)
        
        
   