# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 16:18:59 2023

@author: ZimmermannP
"""

from src.heatex import evap_z
from src.geometrie import geo
from src.air import air
from src.Refrigerant import refrigerant 
import math

R = 'R717' # Kältemittel
ref = refrigerant('R717')

""" Kältemittel/Kreislauf Parameter """
te = -30
M_r = 0.0226 # TODO: Massenstrom oder Leistung als Vorgabe?
#pe  = PropsSI('P','T',273.15+te,'Q',1,R)
lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te,'R717')
x1 = 0 # Start Verdampfung --> kann auch berechnet werden
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
L = 0.1 # Rohrlänge [m]
lR = 220 # Wärmeleitfähigkeit Rippe --> Alu [W/m²K]
lG = 16 # Wärmeleitfähigkeit Grundrohr --> Stahl [W/m²K]
m = M_r / (i_R*i_Q/z*(math.pi/4)*di**2) # Massenstromdichte [kg/m²s]
tg = te # Grundrohrtemperatur tg = te vereinfacht angenommen
""" Erzeugung der Gemotrie und Luft """
geo1 = geo(t_R,sigma_R,di,da,s,sq,sl,L,i_R,i_Q,z,Al,Vl,R_p,lG,lR)
air1 = air(tl, phi_i, Vl)
alpha_i = 680
alpha_a = 31

e1 = evap_z(air1,geo1,alpha_i,alpha_a,tg,m,te)


