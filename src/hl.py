# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 14:18:47 2023

@author: ZimmermannP
"""
import math
#from CoolProp.CoolProp import PropsSI
te = -30
di = 0.0155
rho_g = 1.03
rho_f = 677.73
x = 0.127606
m = 0.0226


m_g = x * m 
print('massenstrom gas', m_g)
V_g = m_g / rho_g
print('Volumenstrom gas',V_g)

m_f =(1-x) * m
print('massenstrom flüssig',m_f)
V_f = m_f / rho_f
print('Volumenstrom flüssig',V_f)

Vg = V_g + V_f
print('V gesamt',Vg)
Ai = math.pi * (di/2)**2 * 6
print('Innenfläche=',Ai)
w = Vg / Ai
print('fließgeschwindigkeit',w)

A1 = V_f / w
h = math.sqrt(2/math.pi) * math.sqrt(A1)
hl = h/di
print(hl)