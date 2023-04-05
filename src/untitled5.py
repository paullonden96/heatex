# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:34:54 2023

@author: ZimmermannP
"""

import random
import math
x = 0.3 
di = 0.014
rho_g = 3.954
rho_l = 633.3
m = 60
g = 9.81
Pi = 3.14
sigma_r  = 0.025391948724494583 
hl = 0.204
phi_wet = math.radians(107)
G1 = 0
G2 = 1
eta_l = 0.000163314462
# Fr_We_l = (g*di**2*rho_l) / sigma_r
Fr_We_l = sigma_r / (g*di**2*rho_l)
a = (x/rho_g) * ((1+0.12*(1-x))*((x/rho_g)+((1-x)/rho_l))+((1.18*(1-x)*\
    (g*sigma_r*(rho_l-rho_g))**0.25)/(m*rho_l**0.5)))**-1

Re_l  = (m * (1-x)*di) / eta_l    
xi = (1.8 * math.log(Re_l) - 1.64)**-2   
Br = 3
print(a)



dis = ((hl*xi*m**2*x**2)/(Br*rho_g*a**(5/2)))**0.5
G_wave = (((g*di*rho_l*rho_g*a**3*Pi)/(2*x**2*math.sqrt(2*(1-math.cos(phi_wet)))))*\
          (1+(Pi**2/(25*hl**2))*(1-x)**G1*(Fr_We_l)**G2))**0.5
    
G_wave_dis = (((g*di*rho_l*rho_g*a**3*Pi)/(2*x**2*math.sqrt(2*(1-math.cos(phi_wet)))))*\
              (1+(Pi**2/(25*hl**2))*(1-x)**G1*(Fr_We_l)**G2+dis))**0.5
    
print('dis=',dis)
print('G_wave=',G_wave_dis)
print('G_wave ohne=',G_wave)