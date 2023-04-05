# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 17:57:24 2023

@author: ZimmermannP

Einfache Verdampfung waagerechtes Rohr 
"""

from src.air import air
from src.geometrie import geo
from src.Refrigerant import refrigerant
from src.model_zürcher import alpha_x
from CoolProp.CoolProp import PropsSI
import math
R = 'R717'
ref = refrigerant('R717')

# tl1 = -18
# tg = -30
# Aa = 114
# a_as = 31
# M_l = 3.84
# dhdt = 1120

# tl2 = tl1 - (tl1-tg) * ( 1- math.exp(-1*((a_as*Aa)/(M_l*dhdt))))
# print(tl2)

""" Luftparameter """
tl = -18 # Lufteintrittstemperatur
phi_i = 0.91 # rel. Luftfeuchte Eintritt
Vl = 10000 #Luftvolumenstrom

t_R = 0.004
sigma_R = 0.00025
di = 0.0155
da = 0.0165
sq = 0.05
sl = 0.05
i_R = 12
i_Q = 4
z = 6
p = i_R * i_Q / z
Al = 1.2
Vl = 10000
R_p = 10*-6

""" Rohrparameter """
tg = -29 # tg = te vereinfacht angenommen
di = 0.0155
da = 0.0165
alpha_i = 680
alpha_a = 31
lw = 16 #Stahl
s = 1/1000
L = 0.1

""" Kältemittelparameter """
te = -30
M_r = 0.0226 # TODO: Massenstrom oder Leistung als Vorgabe?
pe  = PropsSI('P','T',273.15+te,'Q',1,R)
lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te,'R717')


class rohr_z:
    
    def __init__(self,di,da,s,lw,alpha_i,alpha_a,L,tg,M_r):
        self.di = di
        self.da = da 
        self.s = s
        self.lw = lw
        self.alpha_i = alpha_i
        self.alpha_a = alpha_a
        self.L = L
        self.tg = tg
        self.M_r = M_r
        self.m = self.M_r / (p*(math.pi/4)*self.di**2)
        
    
    def k(self,Aa,Ai):
        """Wärmedurchgang"""
        Rg = self.s/self.lw
        ki = 1 / ((Aa/Ai)*(1/self.alpha_i+Rg)+(1/alpha_a))
        return ki
    
    def Q(self,Aa,k):
        """Wärmestrom"""
        global tl
        dT = abs(tl-self.tg)
        Q_a = k*Aa*dT
        return Q_a
    
    def x_z(self,Q):
        """Massendampfgehalt Stelle z"""
        global R
        h1 = PropsSI('H','T',273.15+te,'Q',0,R)
        h2 = PropsSI('H','T',273.15+te,'Q',1,R)
        dh = Q / self.M_r
        x = dh/(h2-h1)
        return x



if __name__ == "__main__":
    luft = air(tl,phi_i,Vl)
    geo1 = geo(t_R,sigma_R,di,da,sq,sl,L,i_R,i_Q,z,Al,Vl,R_p) 
    
    alpha_i = alpha_x
    r0 = rohr_z(di,da,s,lw,alpha_i,alpha_a,L,tg,M_r)
    
    Aa = geo1.A_a()
    Ai = geo1.A_i()
    
    k0 = r0.k(Aa,Ai)
    print(k0)
    Q0 = r0.Q(Aa,k0)
    print(Q0)
    x1 = r0.x_z(Q0)
    print(x1)
    
    
    
    
    

    # Stelle 0 
    
