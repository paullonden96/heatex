# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 18:29:56 2023

Modell für Strömungsform und Wärmeübergang einer 2-Phasen-Strömung beim 
Strömungssieden im horizontalen Rohr nach Dissertation von Olivier Zürcher

-   Strömungskarte stimmt bei te=4 di=14 nicht mit den in der Arbeit abgebildeten
    Diagrammen überein
  
-   Schichten-Wellenströmung konnte nich 100% nachvollzogen werden. Einfluss von 
    Blasensieden auch nicht korrekt.
    
-   Schichtenströmung funkioniert ungefähr +- 20%
     

@author: ZimmermannP
"""

from CoolProp.CoolProp import PropsSI
import math
import matplotlib.pyplot as plt
#from mpmath import coth
import numpy as np
#import pandas as pd
import time
#from sympy import solve,solveset, Symbol, S, sin, Eq
import random
from src.Refrigerant import refrigerant 


Pi = math.pi
g = 9.81

class abschnitt():
    
    def __init__(self,R:str,te,x,m,di,s,theta,lambda_w,q,lambda_g,lambda_l,
                 rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe, dh = False):
        
        self.R = R # Kältemittel
        self.te = te # Verdampfungstemperatur
        self.x = x # Dampfgehalt im Abschnitt
        self.m = m # Massenstromdichte
        self.di = di # Rohrinnendurchmesser
        self.theta = theta # Neigungswinkel zur Horizontalen
        self.q = q # innere Wärmestromdichte
        self.s = s # Wandstärke
        self.lambda_w = lambda_w # Wärmeleitfähigkeit 
        self.lambda_g = lambda_g # Wärmeleitfähigkeit Gas
        self.lambda_l = lambda_l # Wärmeleeitfähigkeit Flüssigkeit
        self.rho_l = rho_l # Dichte Flüssigkeit
        self.rho_g = rho_g # Dichte Gas
        self.eta_g = eta_g # Viskosität Gas
        self.eta_l = eta_l # Viskosität Flüssigkeit
        self.pr_g = pr_g # Prandtl Gas
        self.sigma_r = sigma_r # Oberflächenspannung
        self.cp_g = cp_g # spez. Wärmekapazität Gas
        self.cp_l = cp_l # spez. Wärmekapazität Flüssigkeit
        self.delta_hv = delta_hv # Enthalpiedifferenz Verdampfen 
        self.pe = pe # Verdampfungsdruck
        #self.l = l  # Länge des Abschnittes
        #self.M = M
        if dh == False:
            self.dh = self.di
        else: self.dh = dh
        
        self.Re_l0 = self.m*self.di/self.eta_l
        self.Re_g0 = self.m*self.di/self.eta_g
        self.Re_l = self.m*self.x*self.di/eta_g
        self.Re_g = self.m*(1-self.x)*self.di/self.eta_l
        
        #Martinelli Parameter
        self.X = ((1-self.x)/self.x)**0.875 * (self.rho_g/self.rho_l)**0.5 *\
            (self.eta_l / self.eta_g)**0.125
        #Rouhani 
        self.eps = (self.x/self.rho_g) * ((1+0.12*(1-x))*((x/rho_g)+((1-x)/rho_l))+((1.18*(1-x)*\
            (g*sigma_r*(rho_l-rho_g))**0.25)/(m*rho_l**0.5)))**-1
        
        self.Pr_l = PropsSI('Prandtl','Q',0,'T',te+273.15,R)
        self.Pr_g = PropsSI('Prandtl','Q',1,'T',te+273.15,R)
        

    def h_l_start(self,x):
        ''' Berechnung eines Startwert für Iteration der Flüssigkeitshöhe 
            Einfache Strömungsberechnung unter Annahme homogener Strömung beider Phasen '''
        m_g = x * self.m 
        V_g = m_g / self.rho_g
        m_f =(1-x) * self.m
        V_f = m_f / self.rho_l
        Vg = V_g + V_f
        Ai = math.pi * (self.di/2)**2 * 8
        w = Vg / Ai
        A1 = V_f / w
        h = math.sqrt(2/math.pi) * math.sqrt(A1)
        h_l_start = h/self.di
        if h_l_start >= 1: 
            if self.x < 0.5:h_l_start = 0.6
            if self.x > 0.5:h_l_start = 0.05
        if h_l_start <= 0: h_l_start = 0.001
        #print('h_l_start Wert=',h_l_start)
        return h_l_start

    def strömungsform_x(self,h_l_start):
        h_l_start = h_l_start * self.di
        
        ''' Schichten/Wellenströmung Grenzkurve '''
        a_T, phi_wet = It_a(h_l_start, self.di, self.X) # Dampfvolumenanteil nach Taitel und Duker
        h_l = (self.di/2) * (1- math.cos(phi_wet/2))
        h_l = h_l / self.di
        #print('h_l=',h_l)
        a = a_T
        #print(a)
        cost = 1
        G_strat = (800*((a**2*(1-a))/(self.x**2*(1-self.x)))*self.rho_g*g*(self.rho_l-self.rho_g)*self.eta_l*cost)**(1/3)
        #print('G_strat=',G_strat)
        ''' Wellen/Ringströmung Grenzkurve '''
        G_start = 1000
        a_start = eps(self.x, self.rho_g, self.rho_l, G_start, self.sigma_r)
        Re_l  = (G_start * (1-self.x)*self.di) / self.eta_l
        Fr_We_l = self.sigma_r / (g*self.di**2*self.rho_l)
        q_crit = 0.131 * math.sqrt(self.rho_g) * self.delta_hv * (g*(self.rho_l-self.rho_g)*self.sigma_r)**0.25
        
        if self.q/q_crit >= 0.0374:
            G1 = -24.12*(self.q/q_crit)
            G2 = 4.825*(self.q/q_crit) + 1.053
        else:
            G1 = 0
            G2 = 1
        phi_wet_start, h_l_start = It_p(a_start,self.di)
        G_wavy_start = Wavy(a_start, self.di, self.rho_l, self.rho_g, self.x, phi_wet_start, h_l_start, G1, G2, Fr_We_l, G_start, Re_l)
        G_wavy, a_R = It_Wavy(G_wavy_start, self.x, self.rho_g, self.rho_l, self.sigma_r, self.di, self.eta_l, G1, G2, Fr_We_l,phi_wet_start)
        
        return G_strat, G_wavy, phi_wet, a_T, a_R
    
    def strömungsform(self,G_strat,G_wavy,x_r):
            if self.m <= G_strat:
                f = 'Schichtenströmung'
            elif self.m > G_strat and self.m < G_wavy:
                f = 'Schichten-Wellenströmung'
            elif self.m > G_wavy and self.x < x_r:
                f = 'Intermittet flow'
            elif self.m > G_wavy and self.x > x_r:
                f = 'Ringströmung'
           # TODO: Grenzkurven von Nebel hinzufügen
            else:
                f = None 
            return f
    def alpha_z(self,phi_strat, f,a,G_strat,G_wavy,a_T,alpha_c):
        if f == 'Schichtenströmung':
            alpha_z,phi_wet = a.a_strat(phi_strat,f,a,a_T)
        if f == 'Schichten-Wellenströmung':
            h_l_start = a.h_l_start(self.x)
            alpha_z,phi_wet = a.a_wavy(G_strat,G_wavy,f,h_l_start,a,phi_strat,a_T,alpha_c)
        if f == 'Ringströmung':
            # TODO: Ringströmung hinzufügen
            phi_wet = 2*Pi
            alpha_z = a.a(phi_wet,a,f)
        return alpha_z,phi_wet
                     
    def a_strat(self,phi_wet,f,a,a_T):
        '''local heat transfer stratified flow '''
        C = 0.01361
        m = 0.6965
        phi_s = phi_strat(a_T)
        phi_wet = phi_s
        dh_l, dh_v = dh(f,phi_wet,self.di)
        sigma_c = (self.di/2)*(1-math.cos(phi_wet/2))
        Re_s =  (self.m*(1-self.x)*dh_l) / self.eta_l
        #print('Reynold_s=',Re_s)
        Re_v = (self.m*self.x*dh_v) / self.eta_g
        q_crit = q_0nb_x(self.sigma_r,self.te,self.delta_hv,self.rho_g,self.di,phi_wet,sigma_c,C,m,self.Pr_l,self.lambda_l,Re_s)
        q_l = ql(Re_s,self.lambda_l,dh_l,self.lambda_w,phi_wet,self.di,self.q,self.x,self.m,self.Pr_l,self.eta_l)
        #print('q_l=',q_l)
        #print('q_crit=',q_crit)
        """ Blasensieden """
        if q_l < q_crit:
            hnb = 0
        else:
            #hnb = 0
            hnb = h_nb(self.pe,q_l)
        """ konvektives Sieden """
        #dhlw = self.di*(1-a_T)
        Re_dv = (self.di*self.m)/self.eta_g * (self.x/a_T)
        h_cb_l = 4*C*Re_s**m*self.Pr_l**0.4*(self.lambda_l/dh_l)
        h_cb_v = 0.023*Re_dv**0.8*self.Pr_g**0.4*(self.lambda_g/dh_v)
        # print('sigma_c=', sigma_c)
        # print('q_crit=',q_crit)
        # print('ql=',q_l)
        # print('h_cb_l=',h_cb_l)
        # print('h_cb_v=',h_cb_v)
        # print('h_nb=',hnb)
        
        h_wet = h_w(h_cb_l,hnb)
        alpha_z = h_tp(phi_wet,h_wet,h_cb_v)
        return alpha_z, phi_wet
    def a_wavy(self,G_strat,G_wavy,f,h_l_start,a,phi_wet,a_T,alpha_c):
        ''' local heat transfer stratified-wavy flow '''
        G_sw = (self.m-G_strat) / (G_wavy-G_strat)
        a_strwavy = a.a_sw(h_l_start,G_sw,a_T,G_strat)
        n = 38.127 * self.x**2*(1-self.x)
        phi_wet_dl = G_sw**n # dimensionsloser Winkel phi
        phi_strat = phi_wet
        phi_wet = phi_wet_dl*(-phi_strat)+2*Pi*phi_wet_dl + phi_strat
        #print('phi_wet dim =',phi_wet_dl)
        #print('phi_wet=',phi_wet)
        dh_l, dh_v = dh(f,phi_wet,self.di)
         
        C = 0.01361
        m = 0.6965
        Re_s =  (self.m*(1-self.x)*dh_l) / self.eta_l
        Re_v = (self.m*self.x*dh_v) / self.eta_g
        Re_dv = (self.di*self.m)/self.eta_g * (self.x/a_T)
        sigma_c = sigma_crit(phi_strat,phi_wet,a_strwavy,self.di) 
        q_crit = q_0nb_x(self.sigma_r,self.te,self.delta_hv,self.rho_g,self.di,phi_wet,sigma_c,C,m,self.Pr_l,self.lambda_l,Re_s)
        q_l = ql(Re_s,self.lambda_l,dh_l,self.lambda_w,phi_wet,self.di,self.q,self.x,self.m,self.Pr_l,self.eta_l)
        
        if q_l < q_crit:
            hnb = 0
        else:
            hnb = h_nb(self.pe,q_l)
        
        h_cb_l = 4*C*Re_s**m*self.Pr_l**0.4*(self.lambda_l/dh_l)
        h_cb_v = 0.023*Re_dv**0.8*self.Pr_g**0.4*(self.lambda_g/dh_v)
        
        # print('sigma_c=', sigma_c)
        # print('q_crit=',q_crit)
        # print('ql=',q_l)
        # print('h_cb_l=',h_cb_l)
        # print('h_cb_v=',h_cb_v)
        # print('h_nb=',hnb)
        
        h_wet = h_w(h_cb_l,hnb)
        alpha_z = h_tp(phi_wet,h_wet,h_cb_v)    
        # Dryout ab x = 0.9 laut Zürcher
        if self.x > 0.9:
            #Re_v = (self.m*self.di) / self.eta_g
            h_cb_v = 0.023*Re_dv**0.8*self.Pr_g**0.4*(self.lambda_g/self.di)
            # print('h_cb_v=',h_cb_v)
            # print('alpha_c=',alpha_c)
            a = (h_cb_v - alpha_c) / (0.1**2)
            return a*(self.x-0.9)**2 + alpha_c,phi_wet
        return alpha_z,phi_wet
    def a_a(self,phi_wet,a,f):
        '''local heat transfer annular flow'''
        a_R = eps(self.x, self.rho_g, self.rho_l, self.m, self.sigma_r)
        
        dh_l, dh_v = dh(f,phi_wet,self.di)
        C = 0.01361
        m = 0.6965
        Re_s =  (self.m*(1-self.x)*dh_l) / self.eta_l
        Re_v = (self.m*self.x*dh_v) / self.eta_g
        sigma_c = 2*(self.di/2)*(1-math.sqrt(a_R))
        q_crit = q_0nb_x(self.sigma_r,self.te,self.delta_hv,self.rho_g,self.di,phi_wet,sigma_c,C,m,self.Pr_l)
        
        if self.q < q_crit:
            hnb = 0
        else:
            hnb = h_nb(self.pe,self.q)
        
        h_cb_l = 4*C*Re_s**m*self.Pr_l**0.4*(self.lambda_l/dh_l)
        h_cb_v = 0.023*Re_v**0.8*self.Pr_g**0.4*(self.lambda_g/dh_v)
        
        h_wet = h_w(h_cb_l,hnb)
        alpha_z = h_tp(phi_wet,h_wet,h_cb_v)
        return alpha_z

    def a_sw(self,h_l_start,G_sw,a_T,G_strat):
        """ Dampfvolumenanteil in Schichten-Wellenströmung ist eine Berechnung aus 
        a_Schichten und a_Annular
        einmal nach Rouhani und einmal nach Taitel --> danach in Formel einsetzen """
        
    
        """ Dampfvolumenanteil nach Taitel --> für den Bereich Schichtenströmung """
        """ Eq(3.46) Zürcher Iteration nach h_l Höhe Flüssigkeit"""

        Vv = 1.18*(1-self.x)*((self.sigma_r*g*(self.rho_l-self.rho_g))/self.rho_l**2)**0.25 #weighted mean drift velocity
        C0_Strat = self.rho_l*(((self.x/a_T)-((self.rho_g*Vv)/G_strat))/(self.x*self.rho_l+(1-self.x)*self.rho_g))
        C0_A = 1 + 0.12*(1-self.x)
        C0_wavy = C0_Strat + G_sw*(C0_A-C0_Strat)
        a_sw = (self.x/self.rho_g)*(C0_wavy*((self.x/self.rho_g)+((1-self.x)/self.rho_l))+(Vv/self.m))**-1
        #print('a_sw=',a_sw)
        return a_sw
        
def phi_strat(a):
    phi = 0.0001
    a_start = (2*Pi-phi+math.sin(phi)) / (2*Pi)
    
    while abs(a_start/a) < 0.9999 or abs(a_start/a) > 1.0001:
        phi += 0.00001
        a_start = (2*Pi-phi+math.sin(phi)) / (2*Pi)
        #print(a_start)
    return phi
     
def sigma_crit(phi_strat,phi_wet,a_sw,di):
    """ sigma_crit für stratiffied-wavy  
    gesucht wird Psi und Phi zunächst --> half wetting angle --> 0.5*phi_wet? """    
    Psi = 0.5*phi_wet
    i = 0
    #print('Psi=',Psi)
    Phi_start = 0.5
    #print('Phi_start=',Phi_start)
    a_start = a_geo(Psi,Phi_start)
    #print('a_start=',a_start)
    while abs(a_start/a_sw) < 0.999 or abs(a_start/a_sw) > 1.001:
        i += 1
        if a_start < a_sw:
            Phi_start -= 0.00002
        if a_start > a_sw:
            Phi_start += 0.00002
        a_start = a_geo(Psi,Phi_start)
        #print(abs(a_start/a_sw))
        #print('Phi_neu=',Phi_start)
        #time.sleep(0.001)
        if i == 10**6: quit('keine Lösung Iteration sigma_crit')
    Phi = Phi_start
    
    sigma_c = (di/2)*(1-((math.cos(Phi)+math.cos(Psi))/(1+math.cos(Psi-Phi))))
    #print('sigma_c=',sigma_c)
    return sigma_c

def a_geo(Psi,Phi):
    a = 1 - (1/Pi)*((Psi-0.5*math.sin(2*Psi))-((math.sin(Psi)**2)/(math.sin(Psi-Phi)**2))*\
                    ((Psi-Phi)-0.5*math.sin(2*(Psi-Phi))))   
    return a
           
def a_Taitel(h_l_start,X,rho_g,rho_l,eta_l,x):
    global g
    global Pi
    hl = h_l_start
    A = Pi/4
    Al = 0.25*(Pi-math.acos(2*hl-1)+(2*hl-1)*math.sqrt(1-(2*hl-1)**2))
    Av = 0.25*(math.acos(2*hl-1)-(2*hl-1)*math.sqrt(1-(2*hl-1)**2))
    ul = A / Al 
    uv = A / Av
    Sl = Pi - math.acos(2*hl-1)
    Sv = math.acos(2*hl-1)
    Si = math.sqrt(1-(2*hl-1)**2)
    dl = (4*Av) / Sl
    dv = (4*Av) / (Sl+Si)
    n = 0.2
    m = 0.2
    eq = X**2*((ul*dl)**-n*ul**2*(Sl/Al))-((uv*dv)**-m*uv**2*((Sv/Av)+(Si/Al)+(Si/Av)))
    #print('eq-START=',eq)
    
    while abs(0-eq) > 0.015:
        if eq > 10:
            hl += 0.0001
        if eq > 0.1:
            hl += 0.0000005
        if eq > 0 and eq < 0.1:
            hl += 0.0000005
        if eq < -10:
            hl -= 0.0001
        if eq < -0.1:
            hl -= 0.0000005
        if eq > -0.1 and eq < 0:
            hl -= 0.0000005
        
        eq = eq_Taitel(hl,X)
        #print('eq_IT=',eq)
        time.sleep(0.1)
    #print('hl_ENDE=',hl)
    Av = 0.25*(math.acos(2*hl-1)-(2*hl-1)*math.sqrt(1-(2*hl-1)**2))
    a_T = Av/A
    cost = 1 # cos(theta) winkel zur waagerechten
    
    G_s = (226.3**2*((Al*Av**2*rho_g*(rho_l-rho_g)*eta_l*g*cost)/(Pi**3*x**2*(1-x))))**(1/3)
    return a_T, G_s   
    
def eq_Taitel(hl,X):
    A = Pi/4
    Al = 0.25*(Pi-math.acos(2*hl-1)+(2*hl-1)*math.sqrt(1-(2*hl-1)**2))
    Av = 0.25*(math.acos(2*hl-1)-(2*hl-1)*math.sqrt(1-(2*hl-1)**2))
    ul = A / Al 
    uv = A / Av
    Sl = Pi - math.acos(2*hl-1)
    Sv = math.acos(2*hl-1)
    Si = math.sqrt(1-(2*hl-1)**2)
    dl = (4*Av) / Sl
    dv = (4*Av) / (Sl+Si)
    n = 0.2
    m = 0.2
    return X**2*((ul*dl)**-n*ul**2*(Sl/Al))-((uv*dv)**-m*uv**2*((Sv/Av)+(Si/Al)+(Si/Av)))

def h_tp(phi_wet,h_wet,h_vap):
   return (phi_wet*h_wet+(2*Pi-phi_wet)*h_vap) / (2*Pi)

def h_w(h_cb_l,h_nb):
    return (h_cb_l**3+h_nb**3)**(1/3)

def q_0nb_x(sigma_r,te,delta_hv,rho_g,di,phi_wet,sigma_c,C,m,Pr_l,lambda_l,Re_s):
    '''onset of nuclear boiling at the current vapor quality '''
    r_crit = 0.38*10**-6 # m 
    h_cb = C * Re_s**m * Pr_l**0.4 * (lambda_l/sigma_c)
    q_0nb = (2*sigma_r*(te+273.15)*h_cb) / (r_crit*rho_g*delta_hv)
    #print('q_0nb=',q_0nb)
    return q_0nb

def h_nb(pe,ql):
    """ nuclear boiling Cooper """
    p_c = 113.33 * 10**5
    M0 = 17.03
    return 55 * (pe/p_c)**0.12 * (- math.log(pe/p_c))**-0.55*M0**-0.5*ql**0.67

def ql(Re_s, lambda_l,dh_l,lambda_w,phi_wet,di,q,x,m,Pr_l,eta_l):
    """ durschnittliche Wärmestromdichte flüssig """
    """ 1. """
    """ Biot Zahl """
    if Re_s <= 650:
        Nu_d = 3.66        
    else:
        Nu_d = 4.36 
    C = 0.01361
    m = 0.6965
    Re_s =  (m*(1-x)*dh_l) / eta_l
    alpha_l0 = 4*C*Re_s**m*Pr_l**0.4*(lambda_l/dh_l)
    Bi = (alpha_l0 / lambda_w) * 0.5*(Pi-phi_wet/2)*0.5*di
    # print('Biot=',Bi)
    """ 2. Faktor """
    Fl = 1 / (1+Bi)
    # print('Fl=',Fl)
    p = -0.561*x+0.519
    phi_dry = 2*Pi - phi_wet
    # TODO: q funktioniert noch nicht richtig, keine Ahnung was falsch
    return q / (1-(phi_dry/(2*Pi))*Fl**p)
    #return q / (1-(phi_dry/(2*Pi)))
    
def dh(f,phi_wet,di):
    """ hydraulischer Durchmesser in Abhängigkeit der Strömungsform """
    #print('phi_wet=',phi_wet)
    #print('phi_wet=',round(math.degrees(phi_wet),2))
    if phi_wet < Pi:
        Al = di**2/8 *(phi_wet - math.sin(phi_wet))
        Sl = di*(phi_wet/2)
        Si = di*math.sin(phi_wet/2)
        Av = di**2/8*(2*Pi-2*phi_wet+math.sin(phi_wet))
    else: 
        phi_dry = 2*Pi - phi_wet
        Sl = di*(phi_wet/2)
        Av = (di**2/8)*(2*Pi-2*phi_dry+math.sin(phi_dry))
        Al = ((Pi/4)*di**2) - Av
        Si = di*math.sin(phi_dry/2)
    
    if f == 'Schichtenströmung' or f == 'Schichten-Wellenströmung':
        dh_l = (4*Al) / Sl
        dh_v = (4*Av) / (Sl+Si)
    if f == 'Ringströmung' or f == 'Intermittet flow':
        dh_l = (4*Al)/(Pi*di)
        dh_v = (4*Av) / (Pi*di)
    return dh_l, dh_v
        
def It_a(h_l_start, di, X):
    #Startwerte
    i = 0
    phi_wet_start = 4*math.asin(math.sqrt(h_l_start)/math.sqrt(di))
    Sl = di*(phi_wet_start/2)
    Sv = di*(Pi-phi_wet_start/2)
    Si = di*math.sin(phi_wet_start/2)
    a_start = (2*Pi-phi_wet_start + math.sin(phi_wet_start)) / (2*Pi)
    eq_start = X**2 * (Sl**1.2/((1-a_start)**3)) - ((((Sv+Si)**0.2)/a_start**2)*(((Sv+Si)/a_start)+(Si/(1-a_start))))
    #print(di)
    # print('phi-wet_start=',phi_wet_start)
    # print('Sl=',Sl)
    # print('Sv=',Sv)
    # print('Si=',Si)
    # print('a_start=',a_start)
    # print('eq_start=',eq_start)
    # time.sleep(1)
    while abs(0-eq_start) > 0.002:
        i +=1
        if eq_start <=0:
            phi_wet_start -= 0.000002
        else:
            phi_wet_start += 0.000002
        Sl = di*(phi_wet_start/2)
        Sv = di*(Pi-phi_wet_start/2)
        Si = di*math.sin(phi_wet_start/2)
        a_start = (2*Pi-phi_wet_start + math.sin(phi_wet_start)) / (2*Pi)
        eq_start = X**2 * (Sl**1.2/((1-a_start)**3)) - ((((Sv+Si)**0.2)/a_start**2)*(((Sv+Si)/a_start)+(Si/(1-a_start)))) 
        #time.sleep(0.001)
        #print(eq_start)
        if i == 10**8: quit('keine Iterationslösung It_a gefunden')
    phi_wet_ende = phi_wet_start
    a_ende = (2*Pi-phi_wet_ende + math.sin(phi_wet_ende)) / (2*Pi)
    #print('Dampfvolumenanteil =',a_ende)
    #print('phi_wet = ',phi_wet_ende)
    return a_ende, phi_wet_ende

def eps(x, rho_g, rho_l, m, sigma_r):
    ''' Dampfvolumenanteil nach Rouhani '''
    return (x/rho_g) * ((1+0.12*(1-x))*((x/rho_g)+((1-x)/rho_l))+((1.18*(1-x)*\
        (g*sigma_r*(rho_l-rho_g))**0.25)/(m*rho_l**0.5)))**-1    

def Wavy(a,di,rho_l,rho_g,x,phi_wet,hl,G1,G2,Fr_We_l,G,Re_l):
    ''' Wellengleichung nach Zürcher '''
    if Re_l > 4040:
        xi = (1.8 * math.log10(Re_l) - 1.64)**-2
    if Re_l > 2040 and Re_l <= 3960:
        xi = 5.121552*10**-6*Re_l+2.220595*10**-2
    if Re_l > 1960 and Re_l <= 2040:
        a1 = 1960
        fa1 = 64/960
        f1a1 = -64/1960**2
        b1 = 2040
        fb1 = fa1
        f1b1 = 5.121552 * 10**-6
        Psi1 = (Re_l - b1) * f1a1
        Psi2 = (Re_l - a1 ) * f1b1
        Phi1 = (Re_l - b1)**2 * (3*a1-b1-2*Re_l)
        Phi2 = (Re_l - a1)**2 * (3*b1-a1-2*Re_l)
        xi = fa1/(a1-b1)**3 * (Phi1-Phi2) + (((Re_l-b1)*(Re_l-a1))/(a1-b1)**2) * (Psi1+Psi2)
    if Re_l > 3960 and Re_l <= 4040:
       a1 = 3960
       fa1 = (1.8*math.log10(4040)-1.64)**-2
       f1a1 = 5.121552*10**-6
       b1 = 4040
       fb1 = fa1
       f1b1 = -3.6/(math.log(10)*4040) *(fb1)**(3/2)
       Psi1 = (Re_l - b1) * f1a1
       Psi2 = (Re_l - a1 ) * f1b1
       Phi1 = (Re_l - b1)**2 * (3*a1-b1-2*Re_l)
       Phi2 = (Re_l - a1)**2 * (3*b1-a1-2*Re_l)
       xi = fa1/(a1-b1)**3 * (Phi1-Phi2) + (((Re_l-b1)*(Re_l-a1))/(a1-b1)**2) * (Psi1+Psi2)  
    else:
        xi = 64 / Re_l
    Br = 3
    dis = ((hl*xi*G**2*x**2)/(Br*rho_g*a**(5/2)))
    #print('hl in def=',hl)
    G_wave_dis = (((g*di*rho_l*rho_g*a**3*Pi)/(2*x**2*math.sqrt(2*(1-math.cos(phi_wet)))))*\
              (1+(Pi**2/(25*hl**2))*(1-x)**G1*(Fr_We_l)**G2+dis))**0.5
    
    return G_wave_dis

def It_p(a,di):
    p1 = 0.001
    a1 = (2*Pi-p1 + math.sin(p1)) / (2*Pi)
    #print('a_start=',a1)
    while abs(a1/a) < 0.999 or abs(a1/a) > 1.001:
        p1 += 0.001
        a1 = (2*Pi-p1 + math.sin(p1)) / (2+Pi)
        if p1 > (2*Pi):
            print('keine Iterationslösung bei It_p gefunden')
            exit(0)
        #print(abs(a1/a))
    hl = (di/2)*(1- math.cos(p1/2))
    hl = hl/di
    return p1, hl
    
def It_Wavy(G_wavy_start,x,rho_g,rho_l,sigma_r,di,eta_l,G1,G2,Fr_We_l,phi_wet_start):
    ''' Iteration von der Grenzkurve Wellen '''
    G_vergleich = 0
    G_alt = G_wavy_start
    phi_wet = phi_wet_start
    # print('G_wavy_start=',G_wavy_start)
    i = 0
    j = 0
    while abs(G_alt - G_vergleich) > 1:
        #print(i)
        i += 1
        j += 1
        a = eps(x,rho_g,rho_l,G_alt,sigma_r)
        #dh_l, dh_g = dh('Schichten-Wellenströmung',phi_wet,di)
        #print('G_alt=',G_alt)
        Re_l  = (G_alt * (1-x)*di) / eta_l
        # Re_l = (G_alt*dh_l)/(eta_l)*((1-x)/(1-a))
        #print('Re_l=',Re_l)
        if Re_l <= 650:
            G_neu = (650*eta_l) / ((1-x)*di)
        else:
            phi_wet,hl = It_p(a,di)
            G_neu = Wavy(a, di, rho_l, rho_g, x, phi_wet, hl, G1, G2, Fr_We_l, G_alt, Re_l)
        #print('G_neu=',G_neu)
        if i == 20:
        # manchmal kommt man in eine endlose Iterationsschleife, daher wird nach 5 Durchläufen eine neue zufällige Startzahl gewählt
            G_neu = random.randint(600,1000)
            # print('Random hinzugefügt')
        if i == 40:
        # manchmal kommt man in eine endlose Iterationsschleife, daher wird nach 5 Durchläufen eine neue zufällige Startzahl gewählt
            G_neu = random.randint(200,400)
            # print('Random hinzugefügt')
        if i == 60:
        # manchmal kommt man in eine endlose Iterationsschleife, daher wird nach 5 Durchläufen eine neue zufällige Startzahl gewählt
            G_neu = random.randint(400,600)
            # print('Random hinzugefügt')
        if i == 80:
        # manchmal kommt man in eine endlose Iterationsschleife, daher wird nach 5 Durchläufen eine neue zufällige Startzahl gewählt
            G_neu = random.randint(1200,2000)
            i = 0
            # print('Random hinzugefügt')
            
        G_vergleich = G_alt
        G_alt = G_neu
        #print('Vergleich am Ende',abs(G_vergleich-G_alt))
        if j == 10000: 
            print(f'bei x = {x} keine Lösung bei Iteration G_wavy  gefunden')
            exit(0)
            #break
    G_wavy_ende = G_neu
    #print(f'G_wavy_ende = {G_wavy_ende}')
    return G_wavy_ende, a
def x_A( rho_l, rho_g, eta_l, eta_g):
    ''' Berechnung des Grenzwertes des Massendampfgehaltes ab dem Ringströmung auftreten kann '''
    x1 = 0.001
    x = x1
    dx = dx = (1-x1) / 1000
    X_A = 10**-5
    while abs(X_A/0.34) < 0.99 or abs(X_A/0.34) > 1.001 :
        X_A = ((1-x)/x)**0.875 * (rho_g/rho_l)**0.5 *\
            (eta_l / eta_g)**0.125
        #print('X_A',X_A)
        if X_A > 0.34:
            x += dx
        else:
            x -= dx
    return x

def strömungskarte(di,s,lw,te,m,q,t,R='R717',theta=0):
   
    ''' Strömungskarte nach Zürcher --> verbessert für Ammoniak ausgehend von Steiner VDI
    https://www.researchgate.net/publication/37421459_Evaporation_of_Ammonia_in_a_Smooth_Horizontal_Tube_Heat_Transfer_Measurements_and_Predictions/link/56a2666808aef91c8c0eeafb/download
    '''
    ref = refrigerant(R)
    pe  = PropsSI('P','T',273+te,'Q',1,R) #Verdampfungsdruck
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    x_max = 0.9999
    x1 = 0.01
    dx = (x_max-x1) / t
    val = np.empty(shape=(t,6))
    for i in range (t):
        x = (x1+dx*i)
        print(f'{i}.Durchlauf, Massendampfgehalt= {round(x,4)}')
        a = abschnitt(R,te,x,m,di,s,theta,lw,q,lambda_g,lambda_l,
                       rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe)
        h_l_start = a.h_l_start(x)
        if x !=1: G_strat, G_wavy, phi_wet,a_T,a_R = a.strömungsform_x(h_l_start)
        val[i,0] = x
        val[i,1] = 0
        val[i,2] = G_strat
        val[i,3] = G_wavy
        val[i,4] = a_T
        val[i,5] = a_R
        #print(f'Massenstromdichte Schichten={round(G_strat,2)}, Wellen={round(G_wavy,2)}')
    G_min = min(val[:,3])
    #print('G_min=',G_min)
    x_g = x_A(rho_l,rho_g,eta_l,eta_g)
    a2 = abschnitt(R,te,x_g,m,di,s,theta,lw,q,lambda_g,lambda_l,
                   rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,)
    h_l_start2 = a2.h_l_start(x_g)
    G_strat2, G_wavy2, phi_wet2,a2,a3= a2.strömungsform_x(h_l_start2)
    x = []
    x.append(x_g)
    x.append(G_wavy2)
    return val, x
def alpha_crit(di,s,lw,te,m,R,theta,q,lambda_g,lambda_l,rho_l,rho_g,eta_g,
                  eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe):
    x = 0.9 # kritische Massendampfgehalt ab dem Austrocknung beginnt
    a = abschnitt(R,te,x,m,di,s,theta,lw,q,lambda_g,lambda_l,
                  rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe)
    #Starthöhe Flüssigkeit für Iteration
    h_l_start = a.h_l_start(x)
    #Grenzkurven an der Stelle x
    G_strat, G_wavy, phi_wet,a_T,a1 = a.strömungsform_x(h_l_start)
    #Grenzwert für x damit Ringströmung stattfindet
    x_g = x_A(rho_l,rho_g,eta_l,eta_g)
    #Bestimmung der Strömungsform
    f = a.strömungsform(G_strat, G_wavy, x_g)
    #print(f'Strömungsform: {f}')
    #Bestimmung des lokalen alphas
    alpha_z,phi_wet = a.alpha_z(phi_wet,f,a,G_strat,G_wavy,a_T,1)
    return alpha_z
    
    
    
def alpha_m(di,s,lw,te,m,q,t,R='R717',theta=0):
    """ Wärmeübergangskoeffizient in Abhängigkeit des Massendampfgehaltes x bei festem m
        Unabhängig von der Rohrlänge """
    ref = refrigerant(R)
    pe  = PropsSI('P','T',273.15+te,'Q',1,R) #Verdampfungsdruck
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    x1 = 0.01
    dx = (0.999-x1) / (t-1) # TODO: irgendwas stimmt hier net
    val = np.empty(shape=(t,3))
    alpha_c= alpha_crit(di,s,lw,te,m,R,theta,q,lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,
                         pr_g,sigma_r,cp_g,cp_l,delta_hv,pe)
    for i in range(t):
        x = (x1+dx*i)
        print(f'{i}.Durchlauf --> x = {x}')
        alpha_z,phi_wet = alpha_x(x, di, s, lw, te, m, q, alpha_c)
        val[i,0] = x
        val[i,1] = alpha_z
        val[i,2] = phi_wet
    return val
        
def alpha_x(x,di,s,lw,te,m,q,alpha_c=1,theta=0, R='R717',phi=True):
    """ innerer Wärmeübergang an der Stelle x """
    ref = refrigerant(R)
    pe  = PropsSI('P','T',273.15+te,'Q',1,R) #Verdampfungsdruck
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    
    """ alpha an der Stelle x """
    a = abschnitt(R,te,x,m,di,s,theta,lw,q,lambda_g,lambda_l,
                  rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe)
    #Starthöhe Flüssigkeit für Iteration
    h_l_start = a.h_l_start(x)
    #Grenzkurven an der Stelle x
    G_strat, G_wavy, phi_wet,a_T,a1 = a.strömungsform_x(h_l_start)
    #Grenzwert für x damit Ringströmung stattfindet
    x_g = x_A(rho_l,rho_g,eta_l,eta_g)
    #Bestimmung der Strömungsform
    f = a.strömungsform(G_strat, G_wavy, x_g)

    # print(f'hl_s= {round(h_l_start,4)}')
    # print(f'phi_wet = {round(phi_wet,2)}')
    if x > 0.9 and f == ('Schichten-Wellenströmung' or 'Schichtenströmung'):
        print(f'Strömungsform: {f} Dryout')
    else:
        print(f'Strömungsform: {f}')
    #Bestimmung des lokalen alphas
    alpha_z,phi_wet = a.alpha_z(phi_wet,f,a,G_strat,G_wavy,a_T,alpha_c)
    if phi==True:return alpha_z, math.degrees(phi_wet)
    
    else: return alpha_z
       
       
if __name__ == "__main__":
    di = 0.0155
    s = 1/1000
    lw = 16
    te = -30
    m1 = 15
    m2 = 27
    t = 10
    q = 6500
    x = 0.91
    
    alpha = alpha_x(x,di,s,lw,te,m1,q)
    print(alpha)
    # val1 = alpha_m(di, s, lw, te, m1,q,t)
    # val2 = alpha_m(di, s, lw, te, m2,q,t)
    
    # xd = val1[:,0]
    # a1 = val1[:,1]
    # p1 = val1[:,2]
    # a2 = val2[:,1]
    # p2 = val2[:,2]
    
    # plt.plot(xd, a1, label = ('berechnet m15 nach Zürcher')) # M9 geht weil Schichten
    #plt.plot(xd,a2,label=('berechnet m27 nach Zürcher'))
    # plt.ylim(0)
    # plt.xlim(0,1)
    # plt.xlabel('Massendampfgehalt [-]')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
    # plt.title(f'R717; {te}°C; {q/1000} [kW/m²]')
    # plt.grid(True)
    
    
    """ Grenzwerte für Strömungskarte bei Stelle x """
    #Abschnitt erzeugen
    # a = abschnitt(R,te,x,m,di,s,Ra,theta,lw,q,lambda_g,lambda_l,
    #               rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe)
    #Starthöhe Flüssigkeit für Iteration
    # h_l_start = a.h_l_start(x)
    # Grenzkurven an der Stelle x
    # G_strat, G_wavy, phi_wet,a_T,a1 = a.strömungsform_x(h_l_start)
    #Grenzwert für x damit Ringströmung stattfindet
    # x_g = x_A(rho_l,rho_g,eta_l,eta_g)
    #Bestimmung der Strömungsform
    # f = a.strömungsform(G_strat, G_wavy, x_g)
    # print(f'Strömungsform: {f}')
    
    """  Fehlersuche """
    
    # di = 0.012
    # s = 1/1000
    # te = -20
    # lw = 220
    # x1 = 0.05
    # x2 = 0.25
    # m = 27
    # q = 6500
    
    # a1 = alpha_x(x1, di, s, lw, te, m, q)
    # a2 = alpha_x(x2, di, s, lw, te, m, q)
    
    # print(a1)
    # print(a2)
   
    
    
    