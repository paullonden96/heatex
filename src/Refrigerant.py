# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 17:11:15 2022

@author: ZimmermannP
"""
from __future__ import print_function
from CoolProp.CoolProp import PropsSI
import pandas as pd
import CoolProp
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import SimpleCompressionCycle

class refrigerant():
       
        def __init__(self, name):
             self.name = name 
             
        def stoffwerte(self,te,R):
            lambda_g = PropsSI('L', 'Q',1,'T',273+te,R) #Wärmeleitfähigkeit Gas 
            #print(lambda_g)
            lambda_l = PropsSI('L', 'Q',0,'T',273+te,R) #Wärmeleitfähigkeit flüssig   
            #print(lambda_l)
            rho_l = PropsSI('D', 'Q',0,'T',273+te,R) #Dichte flüssig  [kg/m³]
            #print('rho_l=',rho_l)
            rho_g = PropsSI('D', 'Q',1,'T',273+te,R) #Dichte gas  [kg/m³]
            #print('rho_g=',rho_g)
            eta_g = PropsSI('V','Q',1,'T',273+te,R)  #dynamische Viskosität in [Pa*s]
            eta_l = PropsSI('V','Q',0,'T',273+te,R)  #dynamische Viskosität in [Pa*s]
            pr_g = PropsSI('Prandtl','Q',1,'T',273+te,R) #Prandtl Gas
            sigma_r = PropsSI('I','T',273+te,'Q',0,R) #Oberflächenspannung [N/m]
            #print('sigma=',sigma_r)
            cp_g = PropsSI('C','T',273+te,'Q',1,R) #spez. Wärmekapazität Gas [J/kgK]
            cp_l = PropsSI('C','T',273+te,'Q',0,R) #spez. Wärmekapazität Gas [J/kgK]
            
            
            h1 = PropsSI('H','T',273+te,'Q',0,R) #Siedepunkt
            h2 = PropsSI('H','T',273+te,'Q',1,R) #Taupunkt
            delta_hv =  h2 - h1
            
            return lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv
            
        def dP(self,K, te):
            p1 = PropsSI('P','T',273.15+te,'Q',1,self.name)
            p2 = PropsSI('P','T',273.15+te-K,'Q',1,self.name)
            return abs(p1-p2)
        
        def cycle(self,te,tc,eta_com, dT_sh, dT_sc, Qe):
            
            states = [1,"2s",2,'3"',"4'",4,5,"6'",6,'7"',7]
            states_written = ['Verdichter; Eintritt', 'Verdichter; Austritt,isentrop', 
                              'Verdichter; Austritt', 'Taupunkt Verflüssiger', 
                              'Siedepunkt Verflüssiger', 'Verflüssiger Austritt',
                              'Expansionsventil Eintritt', 'Siedepunkt Verdampfung', 
                              'Verdampfer Eintritt', 'Taupunkt Verdamfer',
                              'Verdampfer Austritt' ]
            
            h = [0,0,0,0,0,0,0,0,0,0,0]
            p = [0,0,0,0,0,0,0,0,0,0,0]
            d = [0,0,0,0,0,0,0,0,0,0,0]
            t = [0,0,0,0,0,0,0,0,0,0,0]
            x = [0,0,0,0,0,0,0,0,0,0,0]
            
            #Verdampfung Ende state 7" 
            h[9] = PropsSI('H','T',273+te,'Q',1, self.name)
            p[9] = PropsSI('P','T',273+te,'Q',1, self.name)
            d[9] = PropsSI('D','T',273+te, 'Q',1,self.name)
            t[9] = te
            
            #Verdampfer Austritt state 7
            h[10] = PropsSI('H','T',273+te+dT_sh,'P',p[9], self.name)
            p[10] = p[9]
            d[10] = PropsSI('D','T',273+te+dT_sh,'P',p[9], self.name)
            t[10] = te+dT_sh
            
            #Verdichtung Eintritt state 1
            h[0] = h[10]
            p[0] = p[9]
            d[0] = d[10]
            t[0] = t[10]
                        
            #Verdichtung Ende, isentrop state 2s
            p3 = PropsSI('P','T',273+tc,'Q',1,self.name) #Kondensationsdruck
            s1 = PropsSI('S','T',273+te+dT_sh,'P',p[9],self.name)
            h2_ideal = PropsSI('H','S',s1,'P',611613.55,'R717')
            h[1] = h2_ideal
            
            p[1] = p3
            d[1] = PropsSI('D','S',s1,'P',p3,self.name)
            t[1] = PropsSI('T','S',s1,'P',p3,self.name)-273
            
            #Verdichtung Ende state 2
            h[2] = h[0] + (h2_ideal - h[0]) / eta_com
            t[2] = PropsSI('T','H',h[2],'P',p3,self.name)-273
            d[2] = PropsSI('D','H',h[2],'P',p3,self.name)
            p[2] = p3
                        
            #Kondensator Eintritt Taupunkt Verflüssiger state 3"
            h[3] = PropsSI('H','T',273+tc,'Q',1,self.name)
            p[3] = p3
            d[3] = PropsSI('D','T',273+tc,'Q',1,self.name)
            t[3] = tc
            
            #Kondensator Siedepunkt state 4'
            h[4] = PropsSI('H','T',273+tc,'Q',0,self.name)
            d[4] = PropsSI('D','T',273+tc,'Q',0,self.name)
            p[4] = p3
            t[4] = tc
            
            #Kondensator Austritt state 4
            h[5] = PropsSI('H','T',273+tc-dT_sc,'P',p3,self.name)
            d[5] = PropsSI('D','T',273+tc-dT_sc,'P',p3,self.name)
            t[5] = tc-dT_sc
            p[5] = p3
            
            #Expansion Eintritt state 5
            h[6] = h[5]
            d[6] = d[5]
            t[6] = t[5]
            p[6] = p3
            
            #Siedepunkt Verdampfung state 6'
            h[7] = PropsSI('H','P',p[9],'Q',0,self.name)
            d[7] = PropsSI('D','T',273+te,'Q',0,self.name)
            t[7] = te
            p[7] = p[9]
            
            #Verdampfer Eintritt state 6
            h[8] = h[5]
            t[8] = te
            p[8] = p[9]
            d[8] = PropsSI('D','H',h[8],'P',p[9],self.name)
            x[8] = PropsSI('Q','H',h[8],'P',p[9],self.name)
         
            
            
            cycle_states = pd.DataFrame({'Prozesspunkt': states_written, 'T [°C]': t, 'p [Pa]':p,
                                         'd [m³/kg]':d, 'h [J/kg]':h, 'x [-]':x
                                         }, index = states )
            return cycle_states
            
        def metrics(self, Qe,te,tc,t1,t2, eta_com):
                p1 = (PropsSI('P','T',273+te,'Q',1, self.name))
                p4 = PropsSI('P','T',273+tc,'Q',1,self.name)
                Qo = PropsSI('H','T',273+te+t1,'P',p1, self.name) - PropsSI('H','T',273+tc-t2,'P',p4,self.name)
                p3 = PropsSI('P','T',273+tc,'Q',1,self.name)
                s1 = PropsSI('S','T',273+te,'Q',1,self.name)
                h2_ideal = PropsSI('H','S',s1,'P',p3,self.name)
                h2 = PropsSI('H','T',273+te+t1,'P',p1, self.name)
                h1 = PropsSI('H','T',273+te,'Q',1, self.name)
                h3 = h1 + (h2_ideal - h1) / eta_com
                P = h3 - h2
                COP = Qo / P
                h6 = PropsSI('H','T',273+tc-t2,'P',p4,self.name)
                q6 = PropsSI('Q','H',h6,'P',p1,self.name)
                m =  Qe / (h2-h6)
                print("x after expansion = " +str(round(q6,3)))
                print("COP = " + str(round(COP,2)))
                print('massflow = ' + str(round(m,3)) +' kg/s with ' + str(Qe/1000) +' kW capacity')
         
        def logph(self, te, tc, dT_sh, dT_sc):
                pp = PropertyPlot(self.name, 'PH', unit_system='EUR',tp_limits='DEF')
                pp.calc_isolines(CoolProp.iT, num=25)
                pp.calc_isolines(CoolProp.iSmass, num=15)
                pp.calc_isolines(CoolProp.iQ, num=11)
                cycle = SimpleCompressionCycle(self.name, 'PH', unit_system='EUR')
                cycle.simple_solve_dt(273.15+te, 273.15+tc, dT_sh, dT_sc, 0.68, SI=True)
                cycle.steps = 50
                sc = cycle.get_state_changes()
                import matplotlib.pyplot as plt
                plt.close(cycle.figure)
                pp.draw_process(sc)
                pp.show()
                
def dpK(dp,te,R='R717'):
    ''' Berechnung des Überhitzungsverlustes '''
    T1 = te+273.15
    pe = PropsSI('P','T',T1,'Q',1,R)
    T2 = PropsSI('T','Q',1,'P',pe-dp,R)
    return T1-T2
                        

if __name__ == "__main__":
   R = 'R717'
   #R = 'R12'
   Ref = refrigerant(R)
   te = -30 #Verdampfungstemperatur
   tc = 10 #Kondensationstemperatur
   dT_sh = 7 #Überhitzung
   dT_sc = 2 #Unterkühlung
   eta_comp = 0.68 #Wirkungsgrad Verdichter
   Q = 27200 #Leistung des Verdampfers 
   m_ref = 0.0226 # Kältemittelmassenstrom in kg/s
   cycle_states = Ref.cycle(te, tc, eta_comp, dT_sh, dT_sc, Q)  
   #Ref.logph(-30,10,7,3)
   print(cycle_states)
   lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l, delta_hv = Ref.stoffwerte(te, R)
   print('lambda_g=',lambda_g)
   print('lambda_l=',lambda_l)
   print('rho_l=',rho_l,'[kg/m³]')
   print('rho_g=',rho_g,'[kg/m³]')
   print('eta_g=',eta_g)
   print('eta_l=',eta_l)
   print('sigma_r=',sigma_r,'[N/m]')
   print('cp_g=',cp_g)
   print('cp_l',cp_l)
   

    
    
    