# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 17:11:15 2022

@author: ZimmermannP
"""
from __future__ import print_function
#import pylab
#import sys
#import CoolProp
#from CoolProp.Plots import SimpleCompressionCycle
#from CoolProp.CoolProp import AbstractState
from CoolProp.CoolProp import PropsSI
#from CoolProp.Plots.Common import BasePlot, PropertyDict, SIunits
from CoolProp.Plots import StateContainer
#import CoolProp.Plots 
#from CoolProp.Plots import Plots
import matplotlib.pyplot as plt
"""
https://tonysyu.github.io/raw_content/matplotlib-style-gallery/gallery.html
"""

R1 = "R744"
R2 = "R404A"
R3 = "R134A"

eta_com = 0.68 # Wirkungsgrad Verdichter
Qe = 27.2 #  kW Leistung Verdampfer
    
class refrigerant:
        
        def __init__(self, name):
             self.name = name
                         
            
        def cycle(self,te,tc,eta_com, t1, t2):
        
            cycle_states = StateContainer()
            
            #Verdampfung Ende State #1
            h1 = PropsSI('H','T',273+te,'Q',1, self.name)
            p1 = (PropsSI('P','T',273+te,'Q',1, self.name))
            d1 = PropsSI('D','T',273+te, 'Q',1,self.name)
            s1 = PropsSI('S','T',273+te,'Q',1,self.name)
            cycle_states[1, 'H'] = h1
            cycle_states[1,'T'] = 273+te
            cycle_states[1,'P'] = p1
            cycle_states[1,'D'] = d1
            cycle_states[1,'S'] = s1
            
            #Verdichtung Eintritt (inkl. Überhitzung) State #2
            h2 = PropsSI('H','T',273+te+t1,'P',p1, self.name)
            p2 = p1
            d2 = PropsSI('D','T',273+te+t1,'P',p1,self.name)
            cycle_states[2, 'H'] = h2
            cycle_states[2,'T'] = 273+te+t1
            cycle_states[2,'P'] = p2
            cycle_states[2,'D'] = d2
            
            # Verdichtung Ende state #3
            p3 = PropsSI('P','T',273+tc,'Q',1,self.name) #Kondensationsdruck
            h2_ideal = PropsSI('H','S',s1,'P',p3,self.name)
            h3 = h1 + (h2_ideal - h1) / eta_com
            
            t3 = PropsSI('T','H',h3,'P',p3,self.name)
            d3 = PropsSI('D','H',h3,'P',p3,self.name)
            cycle_states[3,'H'] = h3
            cycle_states[3,'P'] = p3
            cycle_states[3,'T'] = t3
            cycle_states[3,'D'] = d3
            
            #Kondensator Eintritt Taupunkt Verflüssiger state #4
            h4 = PropsSI('H','T',273+tc,'Q',1,self.name)
            p4 = PropsSI('P','T',273+tc,'Q',1,self.name)
            d4 = PropsSI('D','T',273+tc,'Q',1,self.name)
            cycle_states[4,'H'] = h4
            cycle_states[4,'P'] = p4
            cycle_states[4,'T'] = 273+tc
            cycle_states[4,'D'] = d4
            
            #Kondensator Siedepunkt state#5
            h5 = PropsSI('H','T',273+tc,'Q',0,self.name)
            d5 = PropsSI('D','T',273+tc,'Q',0,self.name)
            
            cycle_states[5,'T'] = 273+tc
            cycle_states[5,'P'] = p4
            cycle_states[5,'D'] = d5
            cycle_states[5,'H'] = h5
            
            #Kondensator Austritt state#6
            h6 = PropsSI('H','T',273+tc-t2,'P',p4,self.name)
            d6 = PropsSI('D','T',273+tc-t2,'P',p4,self.name)
            
            cycle_states[6,'T'] = 273+tc-t2
            cycle_states[6,'P'] = p4
            cycle_states[6,'D'] = d6
            cycle_states[6,'H'] = h6
            
            #Expansion Ende / Verdampfer Eintritt state#7
            d6 = PropsSI('D','H',h6,'P',p1,self.name)
         
            cycle_states[7,'T'] = 273+te
            cycle_states[7,'P'] = p1
            cycle_states[7,'D'] = d6
            cycle_states[7,'H'] = h6
            return print(cycle_states)
    
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
                m =  Qe / ((h2-h6) / 1000)
                print("x after expansion = " +str(round(q6,3)))
                print("COP = " + str(round(COP,2)))
                print('massflow = ' + str(round(m,3)) +' kg/s with ' + str(Qe) +' kW capacity')
         
        def logph(self):   #logp - h Diagramm
                steps = 100
                ps = []
                pt = []
                hs = []
                ht = []
                
                pcrit = PropsSI('pcrit', self.name)
                            
                plow = 0.7
                deltap = (pcrit/100000-plow) / steps
                
                for i in range(steps+1):
                    ps.append(plow+ deltap*i)
                    p = (0.7 + deltap*i)*100000
                    hs.append(PropsSI('H','Q',0,'P',p,self.name) / 1000)
                    
                for i in range(steps+1):
                        pt.append(plow+ deltap*i)
                        p = (0.7 + deltap*i)*100000
                        ht.append(PropsSI('H','Q',1,'P',p,self.name) / 1000)
                
                idx = len(hs) - 1
                hs_rev = []
                ps_rev = []
                
                while (idx >= 0):
                  hs_rev.append(hs[idx])
                  ps_rev.append(ps[idx])
                  idx = idx - 1
                
                p_ges = ps_rev + pt
                h_ges = hs_rev + ht
                
                plt.plot(h_ges, p_ges, label = 'R717')
                plt
                plt.yscale('log')
                plt.ylabel('Druck [bar]')
                plt.xlabel('spez. Enthalpie [kJ/kg]')
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
                plt.grid(True)
                plt.xlim(0,max(h_ges)+150)
                plt.title('log(p)-h Diagramm von R717')
                plt.savefig('logph.png', dpi=300)
                plt.show()
    
    
        def density_over_temperature(self, tmax, tmin):  
            # Veränderung des Dichteverhältnis für verschiedene Temperaturen und Kältemittel
            t01 = 273.15 # 0 °C Verdampfungstemperatur
            t02 = 273.15-10 # °C
            t03 = 273.15-20 # °C
            t04 = 273.15-30 # °C
            t05 = 273.15-40 # °C
            tvar = [t01,t02,t03,t04,t05]
            steps = 100
            tmax = max(tvar)
            tmin = max(tvar) - 50
            dv1 = [] # Dichteverhältnis Ammoniak
            dv2 = [] # R744
            dv3 = [] # R404A
            dv4 = [] # R134A
            
            delta_t = (tmax-tmin) / steps
            t=[]
            for i in range(steps):
                 t.append((tmin-273.15) +(i*delta_t))
            
            
            for j in range(steps):
                    t_neu = (tmin+j*delta_t)
                    dg = PropsSI('D','Q',1,'T',t_neu,self.name)
                    df = PropsSI('D','Q',0,'T',t_neu,self.name)
                    dv1.append(df / dg)
                    dg = PropsSI('D','Q',1,'T',t_neu,R1)
                    df = PropsSI('D','Q',0,'T',t_neu,R1)
                    dv2.append(df / dg)
                    dg = PropsSI('D','Q',1,'T',t_neu,R2)
                    df = PropsSI('D','Q',0,'T',t_neu,R2)
                    dv3.append(df / dg)
                    dg = PropsSI('D','Q',1,'T',t_neu,R3)
                    df = PropsSI('D','Q',0,'T',t_neu,R3)
                    dv4.append(df / dg)
            
            plt.plot(t,dv1, label = 'R717')
            plt.plot(t,dv2, label = 'R744')
            plt.plot(t,dv3, label = 'R404A')
            plt.plot(t,dv4, label = 'R134A')
            plt.ylabel('Dichteverhältnis')
            plt.xlabel('Temperatur in °C')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            plt.grid(True)
            plt.title('Dichteverhältnis in Abhängigkeit der Temperatur verschiedener Kältemittel')
            plt.savefig('Dichteverhältnis.png', dpi=300)
            return plt.show()
        
        
        def dynamic_visosity(self):# dynamische Viskosität in Abhängigkeit der Temperatur nach gas und flüssig getrennt
            steps = 100
            t01 = 273.15 # 0 °C Verdampfungstemperatur
            t02 = 273.15-10 # °C
            t03 = 273.15-20 # °C
            t04 = 273.15-30 # °C
            t05 = 273.15-40 # °C
            tvar = [t01,t02,t03,t04,t05]
            tmax = max(tvar)
            tmin = max(tvar) - 50
            dvf = [] # dynamische viskosität flüssig
            dvg = [] # dynamische viskosität gas
            delta_t = (tmax-tmin) / steps
            t=[]
            for i in range(steps):
                 t.append((tmin-273.15) +(i*delta_t))
            
            for j in range(steps):
                    #h_neu = hmin + (j*delta_h)
                    #t_neu = 0 + (j*delta_q)
                    t_neu = (tmin+j*delta_t)
                    vg = PropsSI('V','Q',1,'T',t_neu,self.name)
                    vf = PropsSI('V','Q',0,'T',t_neu,self.name)
                    dvg.append(vg * 10**(-6))
                    dvf.append(vf * 10**(-6))       
            
            plt.plot(t,dvg, label ="gasförmig")
            plt.plot(t,dvf, label ="flüssig")
            plt.ylabel('Dynamische Viskosität in 10^-6 Pa/s ')
            plt.xlabel('Verdampfungstemperatur °C')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            plt.title('Dynamische Viskosität in Abhängigkeit der Verdampfungstemperatur in °C')
            plt.grid(True)
            plt.savefig('Dynamische_Viskosität.png', dpi=300)
            plt.show()
        
        def tension(self):
            # Abhängigkeit der Oberflächenspannung von der Temperatur
            steps = 100
            t01 = 273.15 # 0 °C Verdampfungstemperatur
            t02 = 273.15-10 # °C
            t03 = 273.15-20 # °C
            t04 = 273.15-30 # °C
            t05 = 273.15-40 # °C
            tvar = [t01,t02,t03,t04,t05]
            tmax = max(tvar)
            tmin = max(tvar) - 50
            ig = [] # Oberflächenspannung gas
            ifl = [] # Oberflächenspannung flüssig
            delta_t = (tmax-tmin) / steps
            t=[]
            for i in range(steps):
                 t.append((tmin-273.15) +(i*delta_t))
            
            for j in range(steps):
                    t_neu = (tmin+j*delta_t)
                    tension_gas = PropsSI('I','Q',1,'T',t_neu,self.name)
                    tension_flüssig = PropsSI('I','Q',0,'T',t_neu,self.name)
                    ig.append(tension_gas)
                    ifl.append(tension_flüssig)
                    
            
            plt.plot(t,ig, label = 'gasförmig')
            plt.plot(t,ifl, label = 'flüssig')
            plt.ylabel('Oberflächenspannung N/m')
            plt.xlabel('Temperatur in °C')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            plt.grid(True)
            plt.title('Oberflächenspannung in Abhängigkeit der Temperatur')
            plt.savefig('Oberflächenspannung.png', dpi=300)
            return plt.show()
        def saturated_steam_volume(self):
        # Abhängigkeit des spezifischen Sattdampfvolumens von der Temperatur
            steps = 100
            t01 = 273.15 # 0 °C Verdampfungstemperatur
            t02 = 273.15-10 # °C
            t03 = 273.15-20 # °C
            t04 = 273.15-30 # °C
            t05 = 273.15-40 # °C
            tvar = [t01,t02,t03,t04,t05]
            tmax = max(tvar)
            tmin = max(tvar) - 50
            vm = []
            delta_t = (tmax-tmin) / steps
            t=[]
            for i in range(steps):
                 t.append((tmin-273.15) +(i*delta_t))
            
            for j in range(steps):
                    #h_neu = hmin + (j*delta_h)
                    #t_neu = 0 + (j*delta_q)
                    t_neu = (tmin+j*delta_t)
                    spez_volumen = PropsSI('D','Q',1,'T',t_neu,self.name)
                    vm.append(1 / spez_volumen)
                
            plt.plot(t,vm, label = self.name)
            plt.ylabel('Spezifisches Sattdampfvolumen m³/kg')
            plt.xlabel('Temperatur in °C')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            plt.title('Spezifisches Sattdampfvolumen in Abhängigkeit der Temperatur')
            plt.savefig('Sattdampfvolumen.png', dpi=300)
            plt.grid(True)
            return plt.show()
        def specific_heat_capacity(self):
            # spezifische Wärmekapazität 
            t01 = 273.15 # 0 °C Verdampfungstemperatur
            t02 = 273.15-10 # °C
            t03 = 273.15-20 # °C
            t04 = 273.15-30 # °C
            t05 = 273.15-40 # °C
            tvar = [t01,t02,t03,t04,t05]
            steps = 100
            tmax = max(tvar)
            tmin = max(tvar) - 50
            cpg = [] # spez. Wk gas
            cpf = [] # spez. Wk flüssig
            delta_t = (tmax-tmin) / steps
            t=[]
            for i in range(steps):
                 t.append((tmin-273.15) +(i*delta_t))
            
            for j in range(steps):
                    t_neu = (tmin+j*delta_t)
                    cg = PropsSI('C','Q',1,'T', t_neu, self.name)
                    cf = PropsSI('C','Q',0,'T',t_neu, self.name)
                    cpg.append(cg / 1000)
                    cpf.append(cf / 1000)    
            
            plt.plot(t,cpg, label = 'gasförmig')
            plt.plot(t,cpf, label = 'flüssig')
            plt.ylabel('Spezifisches Wärmekapazität in kJ / (kgK)')
            plt.xlabel('Temperatur in °C')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            plt.title('Spezifisches Sattdampfvolumen in Abhängigkeit der Temperatur')
            plt.savefig('Wärmekapazität.png', dpi=300)
            plt.grid(True)
            return plt.show()
        
        def thermal_conductivity(self):
        # Wärmeleitfähigkeit
            t01 = 273.15 # 0 °C Verdampfungstemperatur
            t02 = 273.15-10 # °C
            t03 = 273.15-20 # °C
            t04 = 273.15-30 # °C
            t05 = 273.15-40 # °C
            tvar = [t01,t02,t03,t04,t05]
            steps = 100
            tmax = max(tvar)
            tmin = max(tvar) - 50
            wlg = [] # Wärmeleitfähigkeit gas
            wlf = [] # Wärmeleitfähigkeit flüssig
            delta_t = (tmax-tmin) / steps
            t=[]
            for i in range(steps):
                 t.append((tmin-273.15) +(i*delta_t))
            
            for j in range(steps):
                    t_neu = (tmin+j*delta_t)
                    wg = PropsSI('L','Q',1,'T',t_neu,self.name)
                    wf = PropsSI('L','Q',0,'T',t_neu,self.name)
                    wlg.append(wg)
                    wlf.append(wf)
                    
            
            plt.plot(t,wlg, label = 'gasförmig')
            plt.plot(t,wlf, label = 'flüssig')
            plt.ylabel('Wärmeleitfähigkeit W/(mK)')
            plt.xlabel('Temperatur in °C')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            plt.grid(True)
            plt.title('Wärmeleitfähigkeit in Abhängigkeit der Temperatur')
            plt.savefig('Wärmeleitfähigkeit.png', dpi=300)
            plt.show()
        
        
        def deltah(T1, T2, p,R):
            h1 = PropsSI('H','P',p,'T',273+T1,R)
            h2 = PropsSI('H','P',p,'T',273+T2,R)
            delta = abs(h1-h2)
            return delta

