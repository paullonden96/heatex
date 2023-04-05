# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 12:33:06 2023

@author: ZimmermannP
"""
import math
from CoolProp.CoolProp import PropsSI
Pi = math.pi

class geo():
    """ Geometrie Berechnung nach Plank Kap. 12 """
    
    def __init__(self,t_R,sigma_R,di,da,s,sq,sl,L,i_R,i_Q,z,Al,Vl,R_p,lG,lR,theta=0):
        self.t_R = t_R  # Lamellenteilung in m
        self.sigma_R = sigma_R # Lamellendicke in m
        self.di = di  # Innendurchmesser Rohr in m
        self.da = da  # Außendurchmesser Rohr in m
        self.s = s # Rohrwandstärke in m 
        self.sq = sq  # Abstand Rohrachse senkrecht in m
        self.sl = sl  # Abstand Rohrachse waagerecht in m
        self.L = L # Rohrlänge in m
        self.i_R =  i_R  # Anzahl der Rohrreihen
        self.i_Q = i_Q  # Anzahl pro Reihe
        self.Al = Al # Fläche Ventilatorkanal in m²
        self.Vl = Vl # Luftvolumenstrom in m³/h
        self.R_p = R_p # 10**-6 #? technisch glatte Rohre?
        self.lG = lG # Wärmeleitfähigkeit Grundrohr (Stahl)
        self.lR = lR # Wärmeleitfähigkeit Rippen (Alu)
        self.z = z # parallelgeschaltete Rohrstränge
        self.p = self.i_R*self.i_Q / self.z # Anzahl an Rohren durch die KM gleichzeitig läuft = 4*12/6
        self.theta = theta # Winkel zur Horizontalen
    #Grundrohr Wärmeleitwiderstand
    def R_G(self):
        sigma_G = (self.da-self.di)/2*1.1# TODO:muss nach angegebenem Verfahren berechnet werden siehe DKV Blatt 2-24
        #sigma_G nicht korrekt!!! nur zum Test
        return sigma_G / self.lG
    
    #Glattrohroberfläche
    def A_G(self):
        return Pi*self.da*(1 - self.sigma_R / self.t_R) * self.L * self.i_R * self.i_Q
    
    #Rippenoberfläche
    def A_R(self):
        return (2*self.sq*self.sl - (Pi/2)*self.da**2)*(self.L/self.t_R)*self.i_R*self.i_Q
    
    #äußerer Oberfläche
    def A_a(self):
        return (2*self.sq*self.sl + Pi*self.da*(self.t_R-self.sigma_R-(self.da/2)))*\
                (self.L/self.t_R)*self.i_R*self.i_Q
                
    #innere Oberfläche
    def A_i(self):
        return self.di*Pi*self.L*self.i_R*self.i_Q
    
    #Hohlraumanteil    
    def phi(self):
        return 1 - (self.sigma_R/self.t_R) - ((Pi*self.da**2*(self.t_R-self.sigma_R))\
                /(4*self.sq*self.sl*self.t_R))
                  
    #Rohrbündelvolumen in m³
    def V(self):
        return self.L * self.sq * self.sl * self.i_R * self.i_Q
        
    #aquivalenter Durchmesser           
    def d_a_e(self):
        return 4 * self.phi() * self.V() / self.A_a()
    
    #aquivalenter Durchmesser bereift
    # def d_a_e_f():
    #     return ((4*geo.V()*geo.phi_F())/(geo.A_F_G()+geo.A_F_R()))
        
    #mittlere Anströmgeschwindigkeit
    def w_L_m(self):
        wL = self.Vl / (3600*self.Al*self.phi()) 
        return wL 
    
    #scheinbare Wärmeübergangskoeffizient
    def alpha_a_s_h(self,t):
        """alpha_a in der Überhitzungszone"""
        Re_L = (self.w_L_m()*self.d_a_e()) / (PropsSI('V','P',101325,'T',273.15+t,'Air') / PropsSI('D','P',101325,'T',273.16+t,'Air')) # Reynoldszahl
        Pr_L = PropsSI('Prandtl','P',101325,'T',273.15+t,'Air') #Prandtl 
        lambda_L = PropsSI('L','P',101325,'T',273.15+t,'Air')
        alpha_a = 0.31*(lambda_L/self.d_a_e())*Re_L**0.625*Pr_L**(1/3)*(self.d_a_e()/self.sl)**(1/3)
        print(Pr_L)
        print(Re_L)
        print(self.d_a_e()*1000)
        print(self.phi())
        """Rippenwirkungsgrad"""
        rho_R =  1.28 * (self.sq/self.da) *(self.sl/self.sq - 0.2)**0
        h_w = (self.da/2)*(rho_R-1)*(1 + 0.35*math.log(rho_R))
        X = ((2*alpha_a)/(self.sigma_R*self.lR))**0.5*h_w
        eta_R = math.tanh(X) / X
        """scheinbarer Wärmeübergang alpha_as"""
        alpha_as = alpha_a*(self.A_G()/self.A_a()+eta_R*(self.A_R()/self.A_a()))
        return alpha_a,alpha_as,eta_R
    
    def q_i_s(self,Q_h,k_h,A_h,Q):
        A_a_s = self.A_a()-A_h
        return ((Q-Q_h)/A_a_s) * (self.A_a()/self.A_i())
 
if __name__ == '__main__':
    tl = -18
    t_R = 0.004
    sigma_R = 0.00025
    di = 0.0155
    da = 0.0165
    sq = 0.05
    sl = 0.05
    L = 2
    i_R = 4
    i_Q = 12
    z = 6
    Al = 1.2
    Vl = 10000
    R_p = 10*-6
    geo1 = geo(t_R,sigma_R,di,da,sq,sl,L,i_R,i_Q,z,Al,Vl,R_p) 
    
    print('GEOMETRIEPARAMETER:')
    print('')
    print(f'Glattrohroberfläche = {round(geo1.A_G(),2)}[m²]')
    print('Rippenoberfläche =',round(geo1.A_R(),2),'m²')
    print('äußere Oberfläche =',round(geo1.A_a(),2),'m²')
    print('Hohlraumanteil =', round(geo1.phi(),2))
    print('Rohrbündelvolumen =', round(geo1.V(),3), 'm³')
    print('äquivalenter Durchmesser =', round(geo1.d_a_e(),5), 'm')
    print(r'innere Oberfläche=', round(geo1.A_i(),2),'m²')
    
    alpha_a = geo1.alpha_a_s(tl)
    print(f'scheinbarer äußerer Wärmeübergang {round(alpha_a,2)} [W/m²K]')