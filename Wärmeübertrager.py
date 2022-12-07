# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 14:58:28 2022


"""

import math
from CoolProp.CoolProp import PropsSI
Pi = math.pi
"""
äußerer Wärmeübergangskoeffizient Überhittzungszone
ohne Kondensat
Lamellenrohrbündel
flüchtend versetztzt
"""

class heatex():
    Q = 15
    te = -30
    def __init__(self, fluid1, fluid2, ):
        self.fluid1 = fluid1
        self.fluid2 = fluid2
            
verdampfer = heatex('R717','Air')
R = verdampfer.fluid1
RL = verdampfer.fluid2

di = 0.0155 # Innendurchmesser Rohr
da = 0.0165 # Außendurchmesser Rohr
sq = 0.05 # Abstand Rohrachse senkrecht
sl = 0.05 # Abstand Rohrachse waagerecht
L = 2 # Rohrlänge in m
i_R = 32 # Anzahl der Rohrreihen
i_Q = 6 # Anzahl pro Reihe
lambda_stahl = 16 # 
lambda_alu = 220 # W/(mK)
t_R = 0.004 #Lamellenteilung in m
s_R = 0.00025 # Lamellendicke in m
Vl = 10000 # Luftvolumenstrom in m³/h
Al = 1.2 # Fläche Ventilatorkanal in m²
tr = -15 #Lufttemperatur im Raum in °C
te = heatex.te #Verdampfungstemperatur

def deltah(T1, T2, p,R):
    h1 = PropsSI('H','Q',1,'T',273+T1,R)
    h2 = PropsSI('H','P',p,'T',273+T2,R)
    delta = abs(h1-h2)
    return delta


Phi = 1 - (s_R / t_R) - ((Pi*da**2*(t_R - s_R)) / (4*sq*sl*t_R)) #Hohlraumanteil
Aa = (2*sq*sl + Pi*da*(t_R-s_R-(da/2)))*(L/t_R)*i_Q*i_R # Äußere Oberläche für Lamellenrohr
V = (L*sq*sl*i_R*i_Q) # Rohrbündelvolumen
dae = (4*V*Phi)/Aa #Äuvivalenter Durchmesser

"""    
print('Phi = ' + str(round(Phi,4)))
print('Äußere Oberfläche = '+str(round(Aa,3))+ ' m²')       
print('Äquivalenter Durchmesser = ' + str(round(dae,4)) +' m')
print('Rohrbündelvolumen = ' + str(round(V,4))+' m³')
"""

""" 
Durchströmter Kanal weben Lamellenrohrbündel
"""

wL = Vl / (3600*Al) #ungehinderte Luftgeschwindigkeit [m/s]
#print(round(wL,4) , 'freie Anströmgeshwindigkeit in m/s')
wLm = wL / Phi # mittlere Anströmgeschwindigkeiten [m/s]
#print(round(wLm,4), 'mittlere Anströmgeschwindigkeit in m/s')
Re_L = (wLm*dae) / (PropsSI('V','P',101325,'T',273.15+tr,RL) / PropsSI('D','P',101325,'T',273.16+tr,RL)) # Reynoldszahl
#print('Reynoldszahl ' +str(round(Re_L)))
Pr_L = PropsSI('Prandtl','P',101325,'T',273.15-tr,RL) #Prandtl --> ist korrekt
#print('Prandtlzahl ' , round(Pr_L,4))
Nu_D = 0.31*Re_L**0.625*Pr_L**(1/3)*(dae/sl)**(1/3) #Nußelt-Zahl
#print('Nußelt-Zahl' , round(Nu_D,2))

lambda_L = PropsSI('L','P',101325,'T',273.15-tr,RL)
alpha_a = Nu_D * lambda_L / dae
print('äußerer Wärmeübergangskoeffizient' , round(alpha_a,2), ' W/m²K')


"""
Überhitzungszone Kühlmittelseite
7K Überhitzung --> -30 --> -23 °C
"""
m = 0.023 # kg/s 
md = m / (i_Q*(Pi/4)*di**2)
Re_G = md*di / PropsSI('V','Q',1,'T',273+te,R) #turbulente Strömung bei Re > 1000

Pr_G = PropsSI('Prandtl','Q',1,'T',237+te,R)
h = deltah(-30,-23,PropsSI('P','Q',1,'T',273+te,R),R)
lambda_G = PropsSI('L','Q',1,'T',273.15+te,R)
alpha_h_i = 0.0214 * (lambda_G / di) * (Re_G**0.8 - 100) * Pr_G**0.4*(1+(di/L)**(2/3))

print('innerer Wärmeübergangskoeffizient' , round(alpha_h_i,2), ' W/m²K')




