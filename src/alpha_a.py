# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 11:41:17 2023

@author: ZimmermannP

Rechnungsprogramm zur Analytischen Berechnung eines luftbeaufschlagten
Wärmeübertragers bzw. eines geraden Verdampfungsrohres
"""

import math
from src.geometrie import geo
from src.air import air
from src.Refrigerant import refrigerant
from CoolProp.CoolProp import PropsSI


Pi = math.pi
ref = refrigerant('R717')


#Massestrom Kältemittel 
def M(Q,h3,h1):  
    return Q / (h3-h1)

#Wärmestrom Überhitzungszone
def Q_h(M,h3,h2):
    return M * (h3-h2)

#Massenstromdichte
def m_R(M,p,di):
    return M / (p*(Pi/4)*di**2)

#Wärmeübergangsskoeffizient in der Überhitzungszone
def alpha_i_h(m_R,di,te,L):
    Re_g = m_R*di/(PropsSI('V','Q',1,'T',273.15+te,'R717'))
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, 'R717')               
    alpha_i_h = 0.0214 * (lambda_g/di) * (Re_g**0.8 -100) *\
                pr_g**0.4 * (1+(di/L)**(2/3))
    return alpha_i_h                  

def alpha_scheinbar(ag,geo,ar,eta):
    """ scheinbarer Wärmeübergangskoeff. Rippe/Grundrohr """
    alpha = ag * (geo.A_G()/geo.A_a()) + ar*eta*(geo.A_R()/geo.A_a())
    return alpha
#Wärmedurchgang   
def k(alpha_as,alpha_i,geo):
    """ Gl. 12.82 Plank """
    return ((geo.A_a()/geo.A_i())*(1/alpha_i+geo.R_G())+(1/alpha_as))**-1

def t_R_m(t_L_m, t_G_m, eta_R):
    return t_L_m - (t_L_m - t_G_m)*eta_R

def tm(tl1,tr,k,A_a,M,dhdt):
    """ Gl. 12.84 Plank """
    # M = Massenstrom Luft!!
    tl2 = tl1 - (tl1-tr)*(1- math.exp(-1*((k*A_a)/(M*dhdt))))
    return tl2

def tg_m(tl,ks,a_as,tr):
    """ Berechnung der mittleren  Rohrtemp. """
    tg_m = tl - (ks/a_as) * (tl - tr)
    return tg_m

def etaR(geo,alpha_r):
    rho_R =  1.28 * (geo.sq/geo.da) *(geo.sl/geo.sq - 0.2)**0
    h_w = (geo.da/2)*(rho_R-1)*(1 + 0.35*math.log(rho_R))
    X = ((2*alpha_r)/(geo.sigma_R*geo.lR))**0.5*h_w
    eta_R = math.tanh(X) / X
    #alpha_as = alpha_r*(geo.A_G()/geo.A_a()+eta_R*(geo.A_R()/geo.A_a()))
    return eta_R

def alpha_außen(geometrie,air,tg,tl,DEBUG=False):
    """ Berechnung eines äußeren Wärmeübergangs """
    
    alpha_a,alpha_as,eta_R = geometrie.alpha_a_s_h(tl)
    
    alpha_g_s = air.alpha_g_s(tl,tg,alpha_a) # alpha sensibel Grundrohr
    
    t_r_m = t_R_m(tl,tg,eta_R) # mittlere Rippentemp.
    alpha_r_s = air.alpha_g_s(tl,t_r_m,alpha_a)
    eta_R = etaR(geometrie,alpha_r_s)
    
    # 2.Durchlauf mit Rippenwirkungsgrad und neuem alpha
    t_r_m = t_R_m(tl,tg,eta_R)
    alpha_r_s = air.alpha_g_s(tl,t_r_m,alpha_a)
    
    alpha_a = alpha_scheinbar(alpha_g_s, geometrie, alpha_r_s, eta_R)
    if DEBUG == True:
        print(f'alpha_a={alpha_a} \nalpha_as={alpha_as}')
    return alpha_a


if __name__ == "__main__":
    
    R = 'R717' # Kältemittel
    ref = refrigerant('R717')

    """ Kältemittel/Kreislauf Parameter """
    te = -30
    M_r = 0.0226 # TODO: Massenstrom oder Leistung als Vorgabe?
    #pe  = PropsSI('P','T',273.15+te,'Q',1,R)
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
    i_R = 12 # Anzahl Rohrreihen
    i_Q = 4 # Anzahl Rohre pro Rohrreihe
    z = 6 # parallelgeschaltete Rohreihen
    Al = 1.2 # Fläche Ventilatorkanal [m²]
    R_p = 10*-6 # technisch glatte Rohre?
    L = 2 # Rohrlänge [m] --> Teilsegmente für die alpha berechnet werden soll
    lR = 220 # Wärmeleitfähigkeit Rippe --> Alu [W/m²K]
    lG = 16 # Wärmeleitfähigkeit Grundrohr --> Stahl [W/m²K]
    m = 15 # Massenstromdichte [kg/m²s]
    tg = te # Grundrohrtemperatur tg = te vereinfacht angenommen
    """ Erzeugung der Gemotrie und Luft """
    geo1 = geo(t_R,sigma_R,di,da,s,sq,sl,L,i_R,i_Q,z,Al,Vl,R_p,lG,lR)
    air1 = air(tl, phi_i, Vl)

    alpha_i = 680
    #alpha_as = alpha_außen(geo1,air1,tg,tl)
    alpha_as = 80
    
    
    

    
    