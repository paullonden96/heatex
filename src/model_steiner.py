# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:50:09 2023
https://pythonforundergradengineers.com/unicode-characters-in-python.html

Modell nach Steiner VDI
- alpha(Z) in Abhängigkeit der Stelle Z im Rohr

--> heißt: alpha verändert sich bei gleichem m und x je nach dem an welcher Stelle
            Z im Rohr berechnet wird 
- Strömungskarte wurde als eine besser lesbare Form implementiert 
@author: ZimmermannP
"""

from CoolProp.CoolProp import PropsSI
import math
from mpmath import coth
import numpy as np
import time

from src.Refrigerant import refrigerant 


Pi = math.pi
g = 9.81

class abschnitt():
    
    def __init__(self,R:str,te,z,x,m,di,s,Ra,theta,lambda_w,q,lambda_g,lambda_l,
                 rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,l, dh = False):
        
        self.R = R # Kältemittel
        self.te = te # Verdampfungstemperatur
        self.z = z # Stelle des Abschnittes
        self.x = x # Dampfgehalt im Abschnitt
        self.m = m # Massenstromdichte
        self.di = di # Rohrinnendurchmesser
        self.theta = theta # Neigungswinkel zur Horizontalen
        self.q = q # innere Wärmestromdichte
        self.s = s # Wandstärke
        self.lambda_w = lambda_w # Wärmeleitfähigkeit 
        self.Ra = Ra # Rohrrauigkeit
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
        self.l = l  # Länge des Abschnittes
        #self.M = M
        if dh == False:
            self.dh = self.di
        else: self.dh = dh
        
        self.Re_l0 = self.m*self.di/self.eta_l
        self.Re_g0 = self.m*self.di/self.eta_g
        self.Re_g = self.m*self.x*self.di/eta_g
        self.Re_l = self.m*(1-self.x)*self.di/self.eta_l
        
        #Martinelli Parameter
        self.X = ((1-self.x)/self.x)**0.875 * (self.rho_g/self.rho_l)**0.5 *\
            (self.eta_l / self.eta_g)**0.125
        
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
    def grenzkurve_x(self,h_l_start):
        X2_h_start = X2_h(h_l_start)
        h_l = It_X(self.X,X2_h_start,h_l_start)
        #print('hl_ende=',h_l)
        ui = 2*math.sqrt(h_l*(1-h_l))
        if h_l <= 0.5:
            phi_wet = 2*math.asin(ui)
            U_l = phi_wet/2
            U_g = Pi - U_l
            f_l =  (phi_wet - math.sin(phi_wet))/8
            f_g = Pi/4 - f_l
        elif h_l > 0.5:
            phi_dry = 2*math.asin(ui)
            phi_wet = 2*Pi - phi_dry
            U_g = phi_dry/2
            U_l = Pi - U_g
            f_g = (phi_dry - math.sin(phi_dry))/8
            f_l = Pi/4 - f_g
            ui = math.sin(phi_dry/2)
        cosv = 1 #waagerechtes Rohr cos(v)
        
        #Grenzfunktion Schichten/Wellenströmung nach Steiner
        G_strat = (226.3**2*((f_l*f_g**2*self.rho_g*(self.rho_l-self.rho_g)*self.eta_l*g*cosv)/(Pi**3*self.x**2*(1-self.x))))**(1/3)
        
        #Grenzfunktion Wellen/Ring nach Steiner
        Fr_We_l = self.sigma_r/(g*self.di**2*self.rho_l)
        G_wavy = math.sqrt(((16*f_g**3*g*self.di*self.rho_l*self.rho_g)/(Pi**2*self.x**2*math.sqrt(1-(2*h_l-1)**2)))*\
                           ((Pi**2/(25*h_l**2))*Fr_We_l+1))
        
        
        return G_strat,G_wavy,phi_wet
    def x_A(self):
        # Stelle x ab der Ringströmung auftreten kann
        x_A = x_An(self.eta_l,self.eta_g,self.rho_l,self.rho_g)
        return x_A
    
    def strömungsform_x(self,h_l_start,a,x_A):
        G_strat,G_wavy,phi_wet = a.grenzkurve_x(h_l_start)
        if self.m <= G_strat:
            f = 'Schichtenströmung'
        elif self.m > G_strat and self.m < G_wavy:
            f = 'Schichten-Wellenströmung'
        elif self.m > G_wavy and self.x < x_A:
            f = 'Intermittet flow'
        elif self.m > G_wavy and self.x > x_A:
            f = 'Ringströmung'
       # TODO: Grenzkurven von Nebel hinzufügen
        else:
            f = None 
        return f, phi_wet
        
    def q_0nb(self, einlauf=True):
        ''' Wärmestrom am Blasenentstehungspunkt
            wird immer als Einlauf berechnet? siehe Anhang H3.5 VDI '''
        
        Re = self.m*self.dh / self.eta_l  
        #print('Re Blasenbildung=',Re)
        Pr = self.eta_l*self.cp_l / self.lambda_l
        xi = (1.82*math.log10(Re)-1.64)**-2
        z = 10**-5
        d_z = self.di / z
        #print(d_z)
        r_cr = 0.3*10**-6
        Nu_un = ((xi/8)*(Re-1000)*Pr) / (1+12.7*math.sqrt(xi/8)*(Pr**(2/3)-1))   
        if Re < 2300 or Re > 50000: #kein Übergangsbereich 
            if Re < 2300: #laminar
                if d_z >= 1:d_z = 1
                if einlauf is True:
                    Nu_z = 0.455*Pr**(1/3)*math.sqrt(Re*d_z) #q = konst.
                    # Nu_z = 0.332*Pr**(1/3)*math.sqrt(Re*d_z) #tw = konst.
                else:
                     #Nu_z = (4.36**3 + 1.302**3 * Re * Pr * d_z)**(1/3) 
                     Nu_z = (3.66**3+1.077**3*Re*Pr*d_z)**(1/3)            
            elif Re >= 2300: #turbulent
                if einlauf is True:
                    if d_z >= 1: Nu_z = 0.75*Nu_un
                    else: Nu_z = Nu_un*(1+(1/3)*d_z**(2/3))      
                else: Nu_z = Nu_un                   
                    
            alpha_l0 = (Nu_z * self.lambda_l) / self.dh
        if Re > 2300 and Re < 50000: #Übergangsbereich -> max Wert nehmen
            if d_z >= 1: d_z = 1  
            if einlauf is True:
                 Nu_z1 = 0.455*Pr**(1/3)*math.sqrt(Re*d_z) #q=konst.
                 # Nu_z1 = 0.332*Pr**(1/3)*math.sqrt(Re*d_z) #tw = konst.
            else:
                 Nu_z1 = (4.36**3 + 1.302**3 * Re * Pr * d_z)**(1/3) #q=const
                 # Nu_z1 = (3.66**3+1.077**3*Re*Pr*d_z)**(1/3) #tw=const
            if einlauf is True:
                 if d_z >= 1:Nu_z2 = 0.75*Nu_un
                 else:Nu_z2 = Nu_un*(1+(1/3)*d_z**(2/3))
            else:Nu_z2 = Nu_un
            alpha_l0 = (max(Nu_z1, Nu_z2) * self.lambda_l) / self.di
        #print('alpha_l0 Blasensieden=',alpha_l0)   
        return (2*self.sigma_r*(273.15+self.te)*alpha_l0) / (r_cr*self.rho_g*self.delta_hv)

    
    ''' KONVEKTIVES SIEDEN '''
    ''' Schichtstömung und Wellenströmung gesondert beachten!!! '''
    
    ''' Verhältnis lokales alpha(z)_k / alpha_l0 ; Gl. (3) VDI S.942 '''
    
    def alpha_k_l0(self, phi, ausgebildet=False):
            d_h_g = self.dh*((phi-math.sin(phi))\
                             /(phi+2*math.sin(phi/2)))
            epsi = eps(self.x,self.rho_g,self.rho_l,self.m,self.sigma_r)
            #print('eps=',epsi)

            Re_g = (self.m*self.x*d_h_g) / (self.eta_g*epsi)
            #print('Re_g=',Re_g)
            
            Pr_g = (self.eta_g*self.cp_g) / self.lambda_g
            #print('Pr_g=',Pr_g)
            Pr_l = self.eta_l*self.cp_l / self.lambda_l
            #print('Pr_l=',Pr_l)
            
            xi_g0 = (1.82*math.log10(self.Re_g0)-1.64)**-2
            #print(xi_g)
            xi_l0 = (1.82*math.log10(self.Re_l0)-1.64)**-2
            #print(xi_l)
            xi_g = (1.82*math.log10(Re_g)-1.64)**-2
            
            ''' Berechnung der einphasigen Wärmeübergängen alpha_0 ''' 
            def alpha_0(Re,Pr,xi,d,lam): #lam = lambda 
                d_z = self.di / self.z
                if Re < 2300 or Re > 50000:
                    #print('Re außerhalb Übergangsbereich')
                    if Re < 2300:
                        # print('laminare Strömung')
                        if d_z >= 1: d_z = 1
                        # print('d/z=',d_z)
                            #print('hydrodynamisch ausgebildete Strömung')
                        #Nu_z = (4.36**3 + 1.302**3 * Re * Pr * d_z)**(1/3) #q = konst.
                        Nu_z = (3.66**3+1.077**3*Re*Pr*d_z)**(1/3) #tw = konst.
                        if ausgebildet is not True or self.z == 10**-5:
                            # print('here')
                            Nu_z = 0.455*Pr**(1/3)*math.sqrt(Re*d_z) #q = konst.
                            #Nu_z = 0.332*Pr**(1/3)*math.sqrt(Re*d_z) #tw = konst.     
                    elif Re >= 2300:
                        # print('turbulente Strömung') #Bedingung q = konst und tw = konst gleiche Formel
                        Nu_un = ((xi/8)*(Re-1000)*Pr) / (1+12.7*math.sqrt(xi/8)*(Pr**(2/3)-1))
                        Nu_z = Nu_un
                        if ausgebildet is not True or self.z == 10**-5:
                            if d_z >= 1: Nu_z = 0.75*Nu_un   
                            else: Nu_z = Nu_un*(1+(1/3)*d_z**(2/3))
                                
                    return (Nu_z * lam) / d
                if Re > 2300 and Re < 50000:
                    # print('Re im Übergangsbereich --> beide Nu-Zahlen vergleichen')
                    if d_z >= 1: d_z = 1
                    Nu_z1 = (4.36**3 + 1.302**3 * Re * Pr * d_z)**(1/3) # q = konst.
                    # Nu_z1 = (3.66**3+1.077**3*Re*Pr*d_z)**(1/3) #tw = konst.
                    if ausgebildet is not True or self.z == 10**-5:
                        Nu_z1 = 0.455*Pr**(1/3)*math.sqrt(Re*d_z) #q = konst.
                        # Nu_z1 = 0.332*Pr**(1/3)*math.sqrt(Re*d_z) #tw = konst.
                    Nu_un = ((xi/8)*(Re-1000)*Pr) / (1+12.7*math.sqrt(xi/8)*(Pr**(2/3)-1))
                    Nu_z2 = Nu_un
                    if ausgebildet is not True or self.z == 10**-5:
                        if d_z >= 1: Nu_z2 = 0.75*Nu_un    
                        else: Nu_z2 = Nu_un*(1+(1/3)*d_z**(2/3))
                    return (max(Nu_z1, Nu_z2) * lam) / d
    
            # print('Start Berechnung alpha_g0')
            alpha_g0 = alpha_0(self.Re_g0,Pr_g,xi_g0,self.di, self.lambda_g)
            # print('Start Berechnung alpha_l0')
            alpha_l0 = alpha_0(self.Re_l0, Pr_l, xi_l0,self.di, self.lambda_l)
            # print('alpha_l0=',alpha_l0)
            #print('Start Berechnung alpha_g')
            alpha_g = alpha_0(Re_g,Pr_g,xi_g,d_h_g,self.lambda_g)
                        
            alpha_k_z = ((1-self.x)**0.01 * ((1-self.x)+ 1.2*self.x**0.4*\
                        (self.rho_l/self.rho_g)**0.37)**-2.2 +\
                         self.x**0.01 * ((alpha_g0/alpha_l0)*(1+8*(1-self.x)**0.7*\
                        (self.rho_l/self.rho_g)**0.67))**-2)**-0.5
            # print('alpha_k_z=',alpha_k_z)
            return alpha_k_z, alpha_l0, alpha_g    
        
                
    ''' alpha(z)_k Konvektionssieden unter Beachtung des unbenetzten Rohrumfangs '''
    def alpha_z_k(self,s, alpha_lb, phi, alpha_g):
        
        #print('phi_deg=',phi/(2*Pi)*360)
        phi = 2*Pi - phi
        phi = phi/2 # phi = unbenetzter Winkel!!!!! damit eine meistens bessere Benetzung beachtet wird phi/2
        Phi = phi / (2*Pi) #normierte Bogen
        #print('Phi=',Phi)
        a = alpha_g / alpha_lb #Verhältnis der lokalen Wärmeübergangskoeff.
        #print('a=',a)
        M = ((alpha_lb*self.di)/self.lambda_w) *((Pi**2*(self.di+self.s))/(4*self.s))#Rippenkennzahl
        #print('Rippenkennzahl M=',M)
        
        f1 = a*math.sqrt(M)*Phi*(1-Phi)*coth(math.sqrt(M)*(1-Phi))
        f2 = math.sqrt(a*M)*Phi*(1-Phi)*coth(math.sqrt(a*M)*Phi)
        Psi = 1 + Phi*(1-Phi) * ((1-a)**2/a) * (1-((1-(1-a)*Phi)/(f1+f2)))
        #print('Psi=',Psi)
        
        if  M <=1:
            alpha_z_k = alpha_lb*(1-Phi) + alpha_g*Phi
            #print('alpha_z_k=',alpha_z_k)
        else:
            alpha_z_k = (((1-(1-a)*Phi))/Psi) * alpha_lb
            #print('alpha_z_k=',alpha_z_k)
        return alpha_z_k
        
    def alpha_z_b(self,s):    
        ''' alpha(z)_b lokaler Wärmeübergangskoeff. Blasensieden
            nach Steiner VDI '''
        pc = 113.33 #Tabellenwert für Ammoniak
        #pc = 41.6 #Tabellenwert für R12
        #pe = 1.51*100000
        # pe ändern in self.pe!!!!!!
        p_red = (self.pe/100000) / pc
        #print('p* =',p_red)
        d0 = 1*10**-2
        m0 = 100 
        Ra0 = 1*10**-6
        alpha0 = 36640 #Tabellenwert für Ammoniak
        #alpha0 = 3290 #Tabellenwert für R12
        q0 = 150000 #Tabellenwert für Ammoniak
        #q0 = 20000 #Tabellenwert für R12
        M0 = 17.03 #Tabellenwert für Ammoniak
        #M0 = 120.91 #Tabellenwert für R12
        Cf = 0.789 * (M0 / 2.016)**0.11
        if self.R == 'R12': Cf = 1.06
        #print(Cf)
        lambda_ws = self.s*self.lambda_w
        
        #print('lambda_ws= ', lambda_ws)
        p_cr = pc*0.1*100000 # p bei p* = 0.1 in Pa
        
        rho_1 = PropsSI('D', 'Q',0,'P',p_cr,self.R) #Dichte flüssig bei p_cr 
        rho_2 = PropsSI('D', 'Q',1,'P',p_cr,self.R) #Dichte gasförmig bei p_cr 
        sigma_cr = PropsSI('I','P',p_cr,'Q',0,self.R) #Oberflächenspannung
        Pr_cr = PropsSI('Prandtl','P',p_cr,'Q',0,self.R)
        
        q_cr = 0.144 * self.delta_hv * ((rho_1-rho_2)*rho_2)**0.5 *\
            ((g*sigma_cr)/rho_1)**0.25 * Pr_cr**(-0.245)

        if p_red > 0.1:
            q_cr_pb = (3.2 * p_red**0.45 * (1-p_red)**1.2)*q_cr
        elif p_red <= 0.1:
            q_cr_pb = (1.2 * (p_red**0.17 + p_red**0.8))*q_cr
        
        
        if lambda_ws >= 0.7:
            n = 0.9 - 0.36 * p_red**0.13 #für anorganische Stoffe 
            
            alpha_z_b = Cf*(self.q/q0)**n * (2.692*p_red**0.43 +\
                        (1.6*p_red**6.5/(1-p_red**4.4)))*(d0/self.di)**0.5*\
                        (self.Ra/Ra0)**0.133*(self.m/m0)**0.25*\
                        (1-p_red**0.1*(self.q/q_cr_pb)**0.3*self.x)*alpha0
            #return alpha_z_b
        if lambda_ws < 0.7:
            k = 0.675 + 0.325 * math.tanh(3.711*(lambda_ws-3.24*10**-2))
            # print('k=', k)
            if s == 'Schichtenströmung' or 'Wellenströmung':
                psi = 0.46 + 0.4 * math.tanh(3.387*(lambda_ws-8.62*10**-3))
            elif s == 'Schwallströmung':
                psi = 0.671 +  0.329 * math.tanh(3.691*(lambda_ws-8.42*10**-3))
            elif s == 'Ringströmung':
                psi = 0.755 + 0.245 * math.tanh(3,702*(lambda_ws-1.25*10**-2))
            else:
                psi = 1
            n = k*(0.9-0.36*p_red**0.13)
            Cf = psi*Cf
            #print('Cf=',Cf)
            #print('psi=',psi)
            '''Gl.(31) für waagerechte Rohre '''
            alpha_z_b = Cf*(self.q/q0)**n * (2.692*p_red**0.43 +\
                        (1.6*p_red**6.5/(1-p_red**4.4)))*(d0/self.di)**0.5*\
                        (self.Ra/Ra0)**0.133*(self.m/m0)**0.25*\
                        (1-p_red**0.1*(self.q/q_cr_pb)**0.3*self.x)*alpha0
        
        #print('alpha_z_b Gl(31)=',alpha_z_b)
        
        ''' Gl (21) für senkrechte Rohre '''
        n2 = 0.8 - 0.1 * 10**(0.76*p_red)
        Cf2 = 0.435 * (M0 / 2.016)**0.27
        alpha_z_b_senkrecht = Cf2*(self.q/q0)**n2 *\
            (2.816*p_red**0.45 + (3.4 + (1.7/(1-p_red**7)))*p_red**3.7)*\
            (d0/self.di)**0.4 * (self.Ra/Ra0)**0.133 *alpha0
        
        #print('alpha_z_b Gl(21)=',alpha_z_b_senkrecht) 
        if alpha_z_b_senkrecht > alpha_z_b:
            print('Grenze des Einflusses der Massenstromdichte noch nicht erreicht ')
        else:
            print('Massenstromdichte hat keinen Einfluss mehr')
            alpha_z_b = alpha_z_b_senkrecht
        return alpha_z_b
        
    ''' alpha(z) '''
    def alpha_z(self,alpha_z_k, alpha_z_b):
        return (alpha_z_k**3 + alpha_z_b**3)**(1/3)

def It_X(X1,X2,hl):
   # X1 = Martinelli
   # X2 = Martinelli Formel
   X1 = X1**2
   while abs(X2/X1) < 0.999 or abs(X2/X1) > 1.001:
       if X2 > X1:
           hl -= 0.000001
       else:
           hl += 0.000001
       X2 = X2_h(hl)
       #print(abs(X2/X1))
   return hl
                        
            
def alphaz(x,z,di,s,lw,x1,L,R,te,m,Ra,q,theta=0,ausgebildet= False,):
   
    ''' lokales alpha an einer Stelle z berechnen '''
    ref = refrigerant(R)
    pe = PropsSI('P','T',273.15+te,'Q',1,R)      
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)         
     #Bestimmung der notwendigen Wärmestromdichte für den Beginn von Blasensieden
    aq = abschnitt(R,te,z,x,m,di,s,Ra,theta,lw,q,lambda_g,lambda_l,
                   rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,L)   
    q_0nb = aq.q_0nb()
    #print('Wärmestromdichte für Blasenbildung=',q_0nb)
    if q_0nb > q: alpha_z_b = 0  #print('nicht genug Wärmestromdichte für Blasenbildung')
    else: alpha_z_b = 0 # TODO: blasensienden wieder hinzufügen
    
    a = abschnitt(R,te,z,x,m,di,s,Ra,theta,lw,q,lambda_g,lambda_l,
                   rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,L)
    #Bestimmung der notwendigen Wärmestromdichte für den Beginn von Blasensieden
    h_l_start = a.h_l_start(x)
    x_A = a.x_A()
    f,phi = a.strömungsform_x(h_l_start,a,x_A) #Bestimmung der Strömungsform s und des Winkels des unbenetzten Rohrs in rad
    #print('strömungsform=',f)
    #print('benetzter Winkel=',math.degrees(phi))
    
    if f == 'Schichtenströmung' or'Schichten-Wellenströmung': 
            #print('Teilbenetzung beachten')
            alpha_z_k_l0, alpha_l0, alpha_g = a.alpha_k_l0(phi,ausgebildet) #alpha_k und alpha_l und alpha_g an der Stelle 0 im Rohr
            alpha_lb = alpha_z_k_l0 * alpha_l0
            #print('alpha(z)_k / alpha_l0 =', alpha_z_k_l0)
            #print('alpha_l0 =',alpha_l0)
            #print('alpha_lb=',alpha_lb)
            #print('alpha_g=',alpha_g)
            alpha_z_k = a.alpha_z_k(s, alpha_lb, phi, alpha_g)
            alpha_z = alpha_z_k
    else:
            alpha_z_k_l0, alpha_l0, alpha_g = a.alpha_k_l0(phi) #alpha_k und alpha_l und alpha_g an der Stelle 0 im Rohr
            alpha_z_k = alpha_z_k_l0 * alpha_l0
            # print('alpha(z)_k / alpha_l0 =', alpha_z_k_l0)
            # print('alpha_l0 =',alpha_l0)
            # print('alpha_lb=',alpha_lb)
            # print('alpha_g=',alpha_g)
            alpha_z = alpha_lb
   
    return alpha_z
def alpha_m_s(di,s,lw,te,m,t,L,R,Ra,theta,q,pe,y=0,ausgebildet=False): 
    
    ''' Berechnung nach Steiner VDI '''
    
    """ lokales alpha von der Stelle z im Rohr abhängig --> wie alpha(x) Graph? """
    
    ''' Berechnung eines mittleren alphas (alpha1) bei definierter Massenstromdichte
        Rückgabe von alpha und val Matrix Dimension (t,11) --> alle Werte über z
        val[:,0] = x Dampfmassengehalt
        val[:,1] = z Stelle im Rohr
        val[:,2] = alpha_z lokaler Wärmeübergangskoeffizient an der Stelle z
        val[:,3] = alpha_z_b lokaler Wärmeübergangskoeffizient an der Stelle z Blasensieden
        val[:,4] = alpha_z_k lokaler teilbenetzter Wärmeübergangskoeffizient an der Stelle z Konvektionssieden
        val[:,5] = alpha_z_k_l0 Verhältnis alpha_z_k / alpha_l0
        val[:,6] = alpha_l0 einphasiger flüssiger Wärmeübergangskoeffizient bei m = m
        val[:,7] = alpha_g einphasiger gasförmiger Wärmeübergangskoeffizient entsprechend dem Dampfvolumenanteil
        val[:,8] = alpha_lb lokaler Wärmeübergangskoeffizient vollbenetzt Konvektionssieden
        val[:,9] = dp_r Reibungsdruckverlust Stelle z
        val[:,10] = dp_s geodätischer Druckabfall
        val[:,11] = dp_b Beschleunigungsdruckabfall
    '''
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R) 
    pe  = PropsSI('P','T',273.15+te,'Q',1,R)
    x1 = 0.0001
    z1 = 10**-5   
    # if x1 == 0: x1 = 10**-5
    dx = (0.999-x1) / t
    dz = (L-z1) / t
    # leere Matrix für die späteren Diagramme
    val = np.empty(shape=(t+1,12))
    form = []
    
    for i in range(t+1): 
            x = (x1+dx*i) # x an der Stelle z ACHTUNG: lineare Änderung von x vereinfacht angenommen!
            #print(f'{i}.Durchlauf, x = {x}')
            z = (z1+dz*i)
            # print('z=',z)
            
            #Erzeugung eines Rohrabschnittes
            a = abschnitt(R,te,z,x,m,di,s,Ra,theta,lw,q,lambda_g,lambda_l,
                           rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,L)
            #Bestimmung der notwendigen Wärmestromdichte für den Beginn von Blasensieden 
            if i == 0:
                q_0nb = a.q_0nb()
                # print('Wärmestromdichte für Blasenbildung=',q_0nb)
                if q_0nb > q: print('nicht genug Wärmestromdichte für Blasenbildung')
            
            h_l_start = a.h_l_start(x)
            #if x != 1: 
            x_A = a.x_A()
            f, phi = a.strömungsform_x(h_l_start,a,x_A) #Bestimmung der Strömungsform s und des Winkels des unbenetzten Rohrs in rad
            # print(f'Phi_stratiffied = {phi}')
            
            #Entscheidung ob Blasensieden beachtet werden muss
            if q_0nb > q:
                alpha_z_b = 0 
                if f == 'Schichtenströmung' or f == 'Schichten-Wellenströmung': 
                        print('Teilbenetzung beachten')
                        alpha_z_k_l0, alpha_l0, alpha_g = a.alpha_k_l0(phi,ausgebildet) #alpha_k und alpha_l und alpha_g an der Stelle 0 im Rohr
                        alpha_lb = alpha_z_k_l0 * alpha_l0
                        print('alpha(z)_k / alpha_l0 =', alpha_z_k_l0)
                        print('alpha_l0 =',alpha_l0)
                        print('alpha_lb=',alpha_lb)
                        print('alpha_g=',alpha_g)
                        alpha_z_k = a.alpha_z_k(s, alpha_lb, phi, alpha_g)
                        alpha_z = alpha_z_k
                        
                else:
                        alpha_z_k_l0, alpha_l0, alpha_g = a.alpha_k_l0(phi) #alpha_k und alpha_l und alpha_g an der Stelle 0 im Rohr
                        alpha_z_k = alpha_z_k_l0 * alpha_l0
                        alpha_lb = alpha_z_k_l0 * alpha_l0
                        # print('alpha(z)_k / alpha_l0 =', alpha_z_k_l0)
                        # print('alpha_l0 =',alpha_l0)
                        # print('alpha_lb=',alpha_lb)
                        # print('alpha_g=',alpha_g)
                        alpha_z = alpha_lb
                        
            else:
                #print('konvektives Sieden und Blasensieden beachten')
                if f == 'Schichtenströmung' or f == 'Schichten-Wellenströmung':    
                        alpha_z_k_l0, alpha_l0, alpha_g = a.alpha_k_l0(phi) #alpha_k und alpha_l und alpha_g an der Stelle 0 im Rohr
                        alpha_lb = alpha_z_k_l0 * alpha_l0
                        # print('alpha(z)_k / alpha_l0 =', alpha_z_k_l0)
                        # print('alpha_l0 =',alpha_l0)
                        # print('alpha_lb=',alpha_lb)
                        # print('alpha_g=',alpha_g)
                        alpha_z_k = a.alpha_z_k(s, alpha_lb, phi, alpha_g)
                        alpha_z_b = a.alpha_z_b(s)
                        # print('alpha_z_b =',alpha_z_b)
                else:
                        alpha_z_k_l0, alpha_l0, alpha_g = a.alpha_k_l0(phi) #alpha_k und alpha_l und alpha_g an der Stelle 0 im Rohr
                        alpha_z_k = alpha_z_k_l0 * alpha_l0
                        alpha_z_b = a.alpha_z_b(s)
                alpha_z = a.alpha_z(alpha_z_k,alpha_z_b)
                        
            # Befüllen der Matrix Spalten 0 - 8 
            # Zuordnung siehe Variabeln
            # Bsp Spalte 0 = x , Spalte 4 = alpha_z_k
            val[i,0] = x
            val[i,1] = z
            val[i,2] = alpha_z
            val[i,3] = alpha_z_b
            val[i,4] = alpha_z_k
            val[i,5] = alpha_z_k_l0
            val[i,6] = alpha_l0
            val[i,7] = alpha_g 
            val[i,8] = alpha_lb
            form.append(f)
    am = np.mean(val[:,2])    
    return am, val 

def alpha_m2(di,s,lw,te,m,dx,dz,t,x1,z1,L,R,Ra,theta,q,lambda_g,lambda_l,
               rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,M, 
               ausgebildet= False, twisted = False): 
    
    ''' Berechnung eines mittleren alphas (alpha1) bei definierter Massenstromdichte
        Rückgabe von alpha und val Matrix Dimension (t,11)
        val[:,0] = x Dampfmassengehalt nicht linear berechnet!
        val[:,1] = z Stelle im Rohr
        val[:,2] = alpha_z lokaler Wärmeübergangskoeffizient an der Stelle z
        val[:,3] = alpha_z_b lokaler Wärmeübergangskoeffizient an der Stelle z Blasensieden
        val[:,4] = alpha_z_k lokaler teilbenetzter Wärmeübergangskoeffizient an der Stelle z Konvektionssieden
        val[:,5] = alpha_z_k_l0 Verhältnis alpha_z_k / alpha_l0
        val[:,6] = alpha_l0 einphasiger flüssiger Wärmeübergangskoeffizient bei m = m
        val[:,7] = alpha_g einphasiger gasförmiger Wärmeübergangskoeffizient entsprechend dem Dampfvolumenanteil
        val[:,8] = alpha_lb lokaler Wärmeübergangskoeffizient vollbenetzt Konvektionssieden
        val[:,9] = dp_r Reibungsdruckverlust Stelle z
        val[:,10] = dp_s geodätischer Druckabfall
        val[:,11] = dp_b Beschleunigungsdruckabfall
    ''' 
    val = np.empty(shape=(t,3))
    dz = (L-z1) / t
    
    alpham1, val1 = alpha_m_s(di,s,lw,te,m,dx,dz,t,x1,z1,L,R,Ra,theta,q,lambda_g,lambda_l,
                   rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,y=0, 
                   ausgebildet= False, twisted = False)   
    
    Tw = q/alpham1 + (273.15+te) 
    #Tw = -27 + 273.15
    a0 = Pi * di * dz
    h1x = PropsSI('H','T',273.15+te,'Q',x1,R)
    h2 = 0
    for i in range(t+1):
        z = (z1+dz*i)
        
        if i == 0:
            x = x1
            h1 = h1x
        else: 
            h1 = h2
        a = alphaz(x,z,di,s,lw,dx,dz,t,x1,z1,L,R,te,m,Ra,theta,q,lambda_g,lambda_l,
                       rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,y=8, 
                       ausgebildet= False, twisted = False)
        q = a * (Tw - (273.15+te))
        print('q=',q)
        Q = q * a0
        print('Q=',Q)
        h2 = Q/M + h1
        x = PropsSI('Q','H',h2,'P',pe,R)
        if x >= 1 or x < 0: x = 0.999999
        val[i-1,0] = x
        val[i-1,1] = z
        val[i-1,2] = a #alpha_z   
    return val 



            
            
def alpham_m(m1,m2,di,s,lw,dx,dz,t,x1,z1,L,R,te,m,Ra,theta,q,lambda_g,lambda_l,
               rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,y = 0,
               ausgebildet = True, twisted = False):
    ''' Berechnung von alpha_m bei m1 bis m2 '''
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R) 
    mx = []
    ax = []
    for i in range(m1,m2+1,1):
        m = i
        mx.append(i)
        alpha,val = alpha_m2(di,s,lw,dx,dz,t,x1,z1,L,R,te,m,Ra,theta,q,lambda_g,lambda_l,
                       rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,y,ausgebildet, twisted)
        ax.append(alpha)
        print(m,alpha)   
    return mx, ax


def eps(x, rho_g, rho_l, m, sigma_r):
    ''' Gasvolumenanteil nach Rouhani '''
    return (x/rho_g) * ((1+0.12*(1-x))*((x/rho_g)+((1-x)/rho_l))+((1.18*(1-x)*\
        (g*sigma_r*(rho_l-rho_g))**0.25)/(m*rho_l**0.5)))**-1
        

def X2_h(h_l):  
    ''' Martinelli Parameter berechnet aus geometrischen Parametern '''
        
    ui = 2*math.sqrt(h_l*(1-h_l))
    if h_l <= 0.5:
        phi_wet = 2*math.asin(ui)
        U_l = phi_wet/2
        U_g = Pi - U_l
        f_l =  (phi_wet - math.sin(phi_wet))/8
        f_g = Pi/4 - f_l
    elif h_l > 0.5:
        phi_dry = 2*math.asin(ui)
        U_g = phi_dry/2
        U_l = Pi - U_g
        f_g = (phi_dry - math.sin(phi_dry))/8
        f_l = Pi/4 - f_g
        ui = math.sin(phi_dry/2)
        
    return (((U_g+ui)/Pi)**0.25*(Pi**2/(64*f_g**2))*\
                (((U_g+ui)/f_g)+(ui/f_l)))*\
            (Pi/U_l)**0.25 * ((64*f_l**3)/(Pi**2*U_l))            
         
def strömungskarte_s(di,s,lw,te,m,q,t,L,R='R717',Ra=0, theta=0):
    ref = refrigerant(R)
    pe  = PropsSI('P','T',273+te,'Q',1,R) #Verdampfungsdruck
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    x_max = 0.9999
    x1 = 10**-3
    dx = (x_max-x1) / t
    #dz = (L-z1) / t
    x1 = 0.001
    x = 0.2
    z = 1
    val = np.empty(shape=(t,4))
    for i in range(t):
        print(i)
        x = (x1+dx*i)
        print(f'{i}.Durchlauf, Massendampfgehalt= {round(x,4)}')
        a = abschnitt(R,te,z,x,m,di,s,Ra,theta,lw,q,lambda_g,lambda_l,
                   rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,L)
       
        h_l_start = a.h_l_start(x)
        if x !=1: G_strat, G_wavy, phi_wet = a.grenzkurve_x(h_l_start)
        x_A = a.x_A()
        
        val[i,0] = x
        val[i,1] = G_strat
        val[i,2] = G_wavy
        
        
    a2 = abschnitt(R,te,z,x_A,m,di,s,Ra,theta,lw,q,lambda_g,lambda_l,
                   rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,L)
    h_l_start2 = a2.h_l_start(x_A)
    G_strat_A, G_wavy_A, phi_wet_A = a2.grenzkurve_x(h_l_start2)
    
    
    return val, phi_wet,x_A,G_wavy_A
    
def x_An(eta_l,eta_g,rho_l,rho_g):
    X = 0.34
    x = 0.001
    eq = (eta_l/eta_g)**0.125 * ((1-x)/x)**0.875 * (rho_g/rho_l)**0.5

    while abs(eq/X) < 0.999 or abs(eq/X) > 1.01:
        x += 0.001
        #print(x)
        eq = (eta_l/eta_g)**0.125 * ((1-x)/x)**0.875 * (rho_g/rho_l)**0.5
        #print(eq) 
        #time.sleep(0.01)
    return x
     
if __name__ == "__main__":
    """ Parameter """
    params = []
    R = 'R717'
    m = 20
    di = 0.014
    s = 0.93 / 1000
    te = 4
    Ra = 35 * 10**-6
    lw = 26 # Wärmeleitfähigkeit Rohr
    pe  = PropsSI('P','T',273+te,'Q',1,R) #Verdampfungsdruck
    ref = refrigerant(R)
    z = 0.155
    x = 0.1
    t = 10
    L = 6 
    x1 = 0.001
    
    theta = 0 #Winkel zur Horizontalen
    q = 7130  #mittlere Wärmestromdichte
    lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
    params = lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l
    # val, phi_wet,x_A,G_wavy_A = strömungskarte_s(di,s,lw,te,m,t,x1,R,Ra,theta,q,
    #                                           lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe,L)
    a_z = alphaz(x,z,di,s,lw,t,x1,L,R,te,m,Ra,theta,q,*params,delta_hv,pe)
    print('alpha',a_z)
    #print('a_z = ',x_A)
    
    
    