# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 12:30:54 2022

@author: Paulchen
"""

import Refrigerant

#s =  {"Prozesspunkt":[1,"2s",2,'3"',"4'",4,5,"6'",6,'7"',7]}

te = -30 # temperatur evaporation °C
tc = 35 # temperatur condensation °C
eta_com = 0.68 
t1 = 7 #superheat
t2 = 2 #subcooling


Ammoniak = Refrigerant.refrigerant('R717') # Ein neues Kältemittel erstellen 
R = Ammoniak.name
print(R)
Ammoniak.cycle(te,tc,eta_com,t1,t2) # calculat simple cycle with given inputs

Ammoniak.metrics(27.3,te,tc,t1,t2,eta_com) 