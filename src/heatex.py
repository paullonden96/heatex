# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 17:57:24 2023

@author: ZimmermannP

Verdampfungsberechnung
"""

from src.Refrigerant import refrigerant
from src.model_zürcher import alpha_x, alpha_crit
from src.model_steiner import alphaz
from src.alpha_a import alpha_außen
from src.plot import plot_graph

import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import math
import numpy as np


class evap_z():
    
    def __init__(self,air,geo,alpha_i,alpha_a,tg,m,te,R):
        """
        Verbindet alle Konstruktor Klassen miteinander
        Bildet einen Teil eines Verdampfers mit der Geometrie "geo" ab.
        Wird von der Luft "air" umgeben.
        Parameters
        ----------
        air : TYPE
            Luft Objekt.
        geo : TYPE
            Geometrie Objekt.
        alpha_i : TYPE
            innerer Wärmeübergang.
        alpha_a : TYPE
            äußerer Wärmeübergang.
        tg : TYPE
            Grundrohrtemperatur.
        m : TYPE
            Massenstromdichte.
        te : TYPE
            verdampfungstemperatur.

        Returns
        -------
        None.

        """
        self.alpha_i = alpha_i 
        self.alpha_a = alpha_a
        self.tg = tg
        self.m = m 
        self.geo = geo
        self.M_r = (self.geo.i_R*self.geo.i_Q/self.geo.z*(math.pi/4)*self.geo.di**2) * self.m
        self.air = air
        self.te = te
        self.R = 'R717'
    
    
    def k(self):
        """Wärmedurchgang"""
        Rg = self.geo.s/self.geo.lG
        print(Rg)
        k_s = 1 / ((self.geo.A_a()/self.geo.A_i())*(1/self.alpha_i+Rg)+(1/self.alpha_a))
        return k_s
    
    def Q(self,k):
        """Wärmestrom"""
        dT = abs(self.air.tl-self.tg)
        Q_a = k*self.geo.A_a()*dT
        return Q_a
    
    def dx(self,Q):
        """Änderung des Massendampfgehaltes"""
        h1 = PropsSI('H','T',273.15+self.te,'Q',0,self.R)
        h2 = PropsSI('H','T',273.15+self.te,'Q',1,self.R)
        dh = Q / self.M_r
        x = dh/(h2-h1)
        return x
    
def q(Q,Ai):
    qi = Q/Ai
    return qi

def L_evap(geo,air,m,tg,tl,te,R,x1=0,x2=1,f=1,plot=False,comp=False):
    """ Berechnet die Länge die benötigt wird, um von x1 zu x2 zu verdampfen"""
    ref = refrigerant(R)
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te,'R717')
    pe  = PropsSI('P','T',273.15+te,'Q',1,R)
    x_l = []
    L_l = []
    x_l.append(x1)
    L_l.append(0)
    if x1 < 0.01: x1 = 0.01
    Q_g = 0
    L_g = 0
    # alpha_as = 80
    alpha_as = alpha_außen(geo,air,tg,tl)
    # Startwert mit Steiner
    alpha_i_start = f*alphaz(x1,0.01,geo.di,geo.s,geo.lG,x1,1,R,te,m,geo.R_p,0)
    e_start = evap_z(air,geo,alpha_i_start,alpha_as,tg,m,te,R)
    k_start = e_start.k()
    Q_start = e_start.Q(k_start)
    q_start = q(Q_start,geo.A_i())
    i = 1
    # alpha_crit = alpha_i an der Stelle x_crit = 0.9
    alpha_c= f*alpha_crit(geo.di,geo.s,geo.lG,te,m,R,geo.theta,q_start,lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,
                         pr_g,sigma_r,cp_g,cp_l,delta_hv,pe) 
    
    x = x1
    while x <= 1:
        if x > 0.9: f = 1        
        alpha_i = f*alpha_x(x,geo.di,geo.s,geo.lG,te,m,q_start,alpha_c,phi=False)
        e0 = evap_z(air,geo,alpha_i,alpha_as,tg,m,te,R)
        k = e0.k()
        Q = e0.Q(k)
        dx = e0.dx(Q)
        x += dx
        print(f'innere Wärmestromdichte = {q_start}')
        print(f'{i}.Durchlauf \nMassendampfgehalt= {round(x,4)} \nalpha_i= {round(alpha_i,1)}')
        i += 1
        Q_g += Q
        L_g += geo.L
        q_i = q(Q,geo.A_i())
        q_start = q_i

        L_l.append(L_g)
        x_l.append(x)
            
        print(f'Länge des Rohres {round(L_g,2)}[m], Gesamtleistung {round(Q_g/1000,2)}[kW]')
    if plot == True and comp == False:
        title = 'Änderung des Massendampfgehaltes beim Strömungssieden von Ammoniak'
        ylabel = 'Massendampfgehalt [-]'
        xlabel = 'Länge des Rohrs [m]'
        l = (f'Glattrohr Faktor x {f}, {geo.di*1000}mm')
        plot_graph(L_l,x_l,l,ylabel,xlabel,title)
    del x
    if comp == True: return L_l,x_l 
    return L_g
        
def L_compare(geo,air,m,tg,tl,te,R,x1=0,x2=1,f2=1,f3=1,f4=1):
    plot = False
    """ maximal 4 Faktoren miteinander vergleichen 
        Faktor wird an jeder Stelle mit alpha_i multipliziert """
    
    f1 = 1 
    L1,x_1 = L_evap(geo,air,m,tg,tl,te,R,x1,x2,f1,plot,comp=True)
    plt.plot(L1,x_1,label = (f'Glattrohr {geo.di*1000} mit Faktor x {f1}'))
    if f2 != 1:
        L2,x2 = L_evap(geo,air,m,tg,tl,te,R,x1,x2,f2,plot,comp=True)
        plt.plot(L2,x2,label = (f'Glattrohr {geo.di*1000} mit Faktor x {f2}'))
    if f3 != 1:
        L3,x3 = L_evap(geo,air,m,tg,tl,te,R,x1,x2,f3,plot,comp=True)
        plt.plot(L3,x3,label = (f'Glattrohr {geo.di*1000} mit Faktor x {f3}'))
    if f4 != 1:
        L4,x4 = L_evap(geo,air,m,tg,tl,te,R,x1,x2,f4,plot,comp=True)
        plt.plot(L4,x4,label = (f'Glattrohr {geo.di*1000} mit Faktor x {f4}'))

    plt.xlabel('Länge des Rohrs [m]')
    plt.ylabel('Massendampfgehalt x [-]')
    plt.title(f'Strömungssieden Ammoniak horizontales Rohr bei te={te}°C, m={m}kg/m²s')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
    plt.xlim(0)
    plt.ylim(0,1)
    plt.grid(True)
    return None

    

if __name__ == "__main__":
    
    # R = 'R717'
    # ref = refrigerant('R717')

    # """ Luftparameter """
    # tl = -18 # Lufteintrittstemperatur
    # phi_i = 0.91 # rel. Luftfeuchte Eintritt
    # Vl = 10000 #Luftvolumenstrom

    # """ Geometrie Verdampfer """
    # t_R = 0.004
    # sigma_R = 0.00025
    # di = 0.0155
    # da = 0.0165
    # sq = 0.05
    # sl = 0.05
    # i_R = 1
    # i_Q = 1
    # z = 1
    # p = i_R * i_Q / z
    # Al = 1.2
    # Vl = 10000
    # R_p = 10*-6

    # """ Rohrparameter """
    # tg = -30 # tg = te vereinfacht angenommen
    # di = 0.0155
    # da = 0.0165
    # #alpha_i = 680
    # #alpha_a = 31
    # lG = 16 #Stahl
    # lR = 220
    # s = 1/1000
    # L = 0.1

    # """ Kältemittelparameter """
    # te = -30
    # m = 15
    # M_r = (p*(math.pi/4)*di**2) * m # TODO: Massenstrom oder Leistung als Vorgabe?
    # pe  = PropsSI('P','T',273.15+te,'Q',1,R)
    # lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te,'R717')


    # luft = air(tl,phi_i,Vl)
    # geo1 = geo(t_R,sigma_R,di,da,s,sq,sl,L,i_R,i_Q,z,Al,Vl,R_p,lG,lR) 
    # qi = 6500
    # print('qi',qi)
    # x = 0.5
    
    # alpha_as = alpha_außen(geo1,luft,tg,tl)
    # alpha_i = alpha_x(x,di,s,lG,te,m,qi,phi=False)
    # print(alpha_as)
    # print('alpha_i=',alpha_i)
    # e0 = evap_z(luft,geo1,alpha_i,alpha_as,tg,m,te,R)
    # k = e0.k()
    # Q = e0.Q(k)
    # print(Q)
    # dx = e0.dx(Q)
    # print(dx)
    # x1 = 0.01
    # q = 6500
    # alpha_i_start = alphaz(x1,0.01,di,s,lG,x1,1,R,te,m,R_p,q)
    # print(alpha_i_start)
    
    
    
    

    # Stelle 0 
    
    pass