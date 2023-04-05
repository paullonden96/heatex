# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 13:43:04 2023

Konstruktor für Luft
@author: ZimmermannP
"""
from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI
import math

p_atmo = 101325 # Umgebungsdruck in Pa 
class air:
    
    def __init__(self,tl,phi_i,Vl):
        self.R_D = 461.5 # Gaskonstante von Wasserdampf
        self.tl = tl    # Lufttemperatur (Eintritt)
        self.phi_i = phi_i # realtive Luftfeuchte (Eintritt)
        self.Pr_L = PropsSI('Prandtl','P',p_atmo,'T',273.15+tl,'Air') # Prandtl Zahl
        self.delta_h_sub = 2.5 * 10**6 #s pezifische Sublimationsenthalpiedifferenz J/kg
        self.t_L_p = HAPropsSI('D', 'T', 273.15+tl, 'R', phi_i, 'P', p_atmo) - 273.15 # Taupunkt der Luft
        self.lambda_L = PropsSI('L','T',273.15+tl,'P',p_atmo,'Air') # Wärmeleitfähigkeit Luft
        self.beta = 3.46 * 10**-2 # m/s^2 Stoffübergangskoeffizient ANNAHME!!!
        self.D = 22 * 10**-6 # m²/s Diffusuionskoeffizient aus Bsp. S.466
          # D = 0.083 * 1 * ((273.15-v_L)/273.15)**1.81 #Schirmer'sche Gleichung Unterschied von 10^3 zum Bsp - ist was anderes aber hängt damit zusammen
        self.cp_L = PropsSI('C','T',273.15+tl,'P',p_atmo,'Air')
        self.rho_L = PropsSI('D','T',273.15+tl,'P',p_atmo,'Air')
        self.M_L = Vl/3600*self.rho_L
        self.deltap_l = 37.2 #Druckabfall in Pa aus Auslegung
        self.x_L1 = HAPropsSI('W', 'T',273.15+tl,'P',p_atmo,'R',phi_i) #Wassergehalt am Eintritt  
    
    # Fläche Überhitzungszone
    def Aa_h(self,Q_h,te,v_L,dT_sh,k_h):
        #Lufttemp. am Übergang von Überhitzungs- und Verdampfungszone
        t_l_strich = self.tl - Q_h / (self.M_L*self.cp_L)
        # mittlerer Temperaturabstand 
        v_m_h = ((t_l_strich - te) - (v_L-(te+dT_sh)))/\
        math.log((t_l_strich-te)/(t_l_strich-(te+dT_sh)))
        #print(f'Temperaturunterschied = {v_m_h}')
        
        # Fläche Überhitzungszone
        Aa_h = Q_h / (k_h*v_m_h)
        return Aa_h
    
    def alpha_g_s(self,tl,v,alpha_a):
        #Verhältnis gesamter Wärmestrom / sensibler Wärmestrom >= 1
        t_l_strich = tl
        x_g = x(v,1)
        QgQs =  1 + (self.delta_h_sub/self.cp_L) * ((self.x_L1 - x_g)/\
                                                (t_l_strich - v))
        #print(f'Verhältnis QgQs = {QgQs}')
        #gesamter Wärmeübergangskoeffizient r=0.81
        return alpha_a * QgQs**0.81
    
    def dhdt(self,Q_h,geo,tR,tG):
        t_l_strich = self.tl - Q_h / (self.M_L*self.cp_L)
        xR = x(tR,1)
        xG = x(tG,1)
        x1 = self.x_L1
        arag = geo.A_a() / geo.A_G()
        
        dxdt = ((x1- xG) / (t_l_strich-tG)) * ((1+arag*((x1-xR)/(x1-xG)))/\
                                               (1+arag*((t_l_strich-tR)/(t_l_strich-tG))))
        dhdt = self.cp_L + dxdt*self.delta_h_sub
        return dhdt
    
#Wassergehalt Luft in kg/kg
def x(t,phi):
     return HAPropsSI('W', 'T',273.15+t,'P',101325,'R',phi)   
            
    
    # #Entfeuchtung pro Kelvin in kg/kgK gute Darstellung S.456
    # def entf(self,t_G, t_R, t_L,geo):
    #     xl = x(t_L, self.phi_i) #Wassergehalt Luft Eintritt
    #     xg = x(t_G, 1) #Wassergehalt Grundrohrtemperatur
    #     xr = x(t_R,1) #Wassergehalt Rippentemperatur
    #     return ((xl-xg)/(t_L-t_G)) * (1+(geo.A_R()/geo.A_G())*((xl-xr)/(xl-xg)))/\
    #            (1+(geo.A_R()/geo.A_G())*((t_L-t_R)/(t_L-t_G)))
    
    # #Reiffdicke in Abhängigkeit der Zeit t in h und v in °C
    # def sigma_f(t,v):
    #     lambda_eis = 2.23 #W/m*K Wärmeleitfähigkeit Eis
    #     F1 = 0.39 #Alternativ 0.435 siehe S.293
    #     rho_eis = 917 #kg/m³
    
    #     return F1 * math.sqrt(lambda_eis / (air.delta_h_sub * rho_eis ))* \
    #            math.sqrt(t*3600 * (air.t_F_o(v) - v))* \
    #             ((air.p_D_L() - air.p_D_F_strich(v)) / (air.p_D_L_strich() - air.p_D_F_strich(v)))**0.27 
   
           
    # #mittlere Dampfteildruck Luft              
    # def p_D(v):
    #     return  4.689 * (1.486 + v/ 100)**12.3 
    
    # #Temperatur an der Reifoberfläche 
    # def t_F_o(v):
    #     return air.t_L_p - 0.22 * (air.t_L_p - v) #return in °C
    
    # # Wasserdampfsättigungsdruck in der Luft in Pa
    # def p_D_L_strich():
    #     return 4.689 * (1.486 + air.v_L / 100)**12.3 
    
    # # Wasserdampfpartialdruck in der Luft in Pa
    # def p_D_L():
    #     return air.phi_i * air.p_D_L_strich() 
    
    # # Wasserdampfsättigungsdruck an der Reifoberfläche in Pa
    # def p_D_F_strich(v):
    #     return 4.689 * (1.486 + (air.t_F_o(v)) / 100)**12.3 
        
    # #mittlere Dampfteildruck im inneren des Reifs
    # def p_D_F_quer(v):
    #      return air.p_D_strich(v) + 0.65 *(air.p_D_F_strich(v) - air.p_D_strich(v))
    
    # #Sättigungsdruck Luft bei Temp. v
    # def p_D_strich(v):
    #     return 4.689 * (1.486 + v/ 100)**12.3
    
    # #Sättigungstemperatur bei p
    # def t_f(p):
    #     return (-61.125785 + 8.1386*(ln(p)) - ((7.422003 *10**(-2)) * \
    #     (ln(p))**2) + (6.283721 *10**(-2)) *\
    #     ((ln(p))**3) - (2.7237063*10**-3 * (ln(p))**4))
        
    # #mittlere Temperatur Luft (Sättigungstemperatur v_s nach Glück / Lufttemp) 
    # def t_quer(v):
    #     v_s = air.t_f(air.p_D_F_quer(v))
    #     return (v_s + air.v_L) / 2
    
    # #delta h / delta t Luft Enthalpieänderung
    # def dhdt(dxdt):
    #     return air.cp_L + dxdt * air.delta_h_sub
    
    # #delta x / delta t Luft
    # def dxdt(t_g,t_l_strich,t_r):
    #     arag = geo.A_R() / geo.A_G()
    #     return ((air.x_L1- air.x(t_g,1)) / (t_l_strich-t_g)) * \
    #         ((1+arag*((air.x_L1-air.x(t_r,1))/(air.x_L1-air.x(t_r,1))))/\
    #          1 + arag * ((t_l_strich-t_r)/(t_l_strich-t_g)))
    
    # # t_L_m mittlere Lufttemperatur
    # def t_L_m(k, t_R_m, dhdt, t_L_m, t_G_m):
    #     return air.v_L - (air.v_L - t_R_m) * ( 1 - math.exp((k*geo.A_a())/(2*air.M_L*air.dhdt())))
    

        
        
    