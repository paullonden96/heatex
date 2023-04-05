# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 12:13:24 2023

@author: ZimmermannP
"""
from src.Refrigerant import refrigerant 
from src.model_zürcher import alpha_m

from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

R = 'R717'
m = 5
di = 0.0155
s = 0.5 / 1000
te = -30
Ra = 35 * 10**-6
lw = 16 # Wärmeleitfähigkeit Rohr
pe  = PropsSI('P','T',273+te,'Q',1,R) #Verdampfungsdruck
ref = refrigerant(R)
x = 0.312
theta = 0 #Winkel zur Horizontalen
q = 5400  #mittlere Wärmestromdichte
t = 20

di = 0.0155
s = 0.5 / 10**3
m = 5
lw = 16
te = -30
q = 5400
t = 15
lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)

val = alpha_m(di, s, lw, te, m, t, R, Ra, theta, q, lambda_g, lambda_l, rho_l, rho_g, eta_g, eta_l, pr_g, sigma_r, cp_g, cp_l, delta_hv, pe)

xd = val[:,0]
alpha_z = val[:,1]

plt.plot(val[:,0], val[:,1], label = (f'Stahl {lw} [W/mK], {m} [kg/m²s] {di*1000:.2f} x {s*1000} [mm] Glattrohr'))
plt.ylim(0,max(alpha_z)+100)
plt.xlim(0,1)
plt.ylabel('\u03B1_i [W/m²K]')
plt.xlabel('Massendampfgehalt [-]')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
#plt.title(f'{R}; {te}°C;')
plt.title(f'{R};  te={te}°C; q={q/1000} kW/m² Glattrohr')
plt.grid(True)