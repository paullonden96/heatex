# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 16:06:26 2023

@author: ZimmermannP
"""
""" Test Iteration Grenzkurve """
m = 1000
a = eps(x, rho_g, rho_l, m, sigma_r)
print('a=',a)
Re_l = (m * (1-x)*di) / eta_l
print('Re_l=',Re_l)
q_crit = 0.131 * math.sqrt(rho_g) * delta_hv * (g*(rho_l-rho_g)*sigma_r)**0.25
print('qcrit=',q_crit)
#print(20000/q_crit)
G1 = 0
G2 = 1
Br = 3
phi_wet,hl = It_p(a,di)
print('hl=',hl)
print('phi_wet=',phi_wet)
xi = (1.8 * math.log(Re_l) - 1.64)**-2
Fr_We_l = sigma_r / (g*di**2*rho_l)  
print('Fr_We_l',Fr_We_l)
dis = ((hl*xi*m**2*x**2)/(Br*rho_g*a**(5/2)))
m1 = math.sqrt(((g*di*rho_l*rho_g*a**3*Pi)/(2*x**2*math.sqrt(2*(1-math.cos(phi_wet)))))*\
            (1+(Pi**2/(25*hl**2))*(1-x)**G1*(Fr_We_l)**G2+dis))
print('m=',m1)

m2 = math.sqrt(((g*di*rho_l*rho_g*a**3*Pi*math.sin(phi_wet/2))/(2*x**2*(1-math.cos(phi_wet))))*\
               (1+(Pi**2/(25*hl**2))*(1-x)**G1*(Fr_We_l)**G2+dis))
print('m2=',m2)

f_g = a*(Pi*(di/2)**2) / di**2
m3 = math.sqrt(((16*f_g**3*g*di*rho_l*rho_g)/(Pi**2*x**2*math.sqrt(1-(2*hl-1)**2)))*\
               (1+(Pi**2/(25*hl**2))*(1-x)**G1*(Fr_We_l)**G2+dis))
    
print('m3=',m3)
