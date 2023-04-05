# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 12:27:18 2023

@author: ZimmermannP
"""

from src.heatex import  L_evap,L_compare
from src.air import air
from src.geometrie import geo
from src.Refrigerant import refrigerant

import math


R = 'R717' # Kältemittel
ref = refrigerant('R717')

""" Kältemittel/Kreislauf Parameter """
te = -30
# M_r = 0.0226 # TODO: Massenstrom oder Leistung als Vorgabe?
# pe  = PropsSI('P','T',273.15+te,'Q',1,R)
lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te,'R717')
x1 = 0.127 # Start Verdampfung --> kann auch berechnet werden
x2 = 1
""" Luftparameter """
tl = -18 # Lufteintrittstemperatur [°C]
phi_i = 0.91 # rel. Luftfeuchte Eintritt [-]
Vl = 10000 # Luftvolumenstrom [m³/h]

""" Geometrie Verdampfer """
t_R = 0.004 # Lamellenteilung [m]
sigma_R = 0.00025 # Lamellendicke [m]
di = 0.0155 # Innendurchmesser Grundrohr [m]
da = 0.0165 # Außendurchmesser Grundrohr [m]
s = 0.001 # Wandstärke Grundrohr [m]
sq = 0.05 # Abstand Rohrachsen Verdampferrohre senkrecht [m]
sl = 0.05 # Abstand Rohrachsen Verdampferrohre waagerecht [m]
i_R = 1 # Anzahl Rohrreihen
i_Q = 1 # Anzahl Rohre pro Rohrreihe
z = 1 # parallelgeschaltete Rohreihen
Al = 1.2 # Fläche Ventilatorkanal [m²]
R_p = 10*-6 # technisch glatte Rohre?
L = 0.2# Rohrlänge [m] --> Teilsegmente für die alpha berechnet werden soll
lR = 220 # Wärmeleitfähigkeit Rippe --> Alu [W/m²K]
lG = 16 # Wärmeleitfähigkeit Grundrohr --> Stahl [W/m²K]
m = 10 # Massenstromdichte [kg/m²s]
tg = te # Grundrohrtemperatur tg = te vereinfacht angenommen
""" Erzeugung der Gemotrie und Luft """
geo1 = geo(t_R,sigma_R,di,da,s,sq,sl,L,i_R,i_Q,z,Al,Vl,R_p,lG,lR)
air1 = air(tl, phi_i, Vl)

""" Rechnung wie lang L sein muss für Verdampfung von x1 - x2 Glattrohr """
#L_gesamt = L_evap(geo1,air1,m,tg,tl,te,R,x1,x2,f=1.5,plot=True)
#print(f'Gesamtlänge= {round(L_gesamt,2)}[m]')

""" Vergleich Strömungssieden mit verschiedenen Multiplikatoren (von x=0-0.9) für alpha_i """
f = 1.5
L_compare(geo1, air1, m, tg, tl, te, R,x1,x2,f,2)


