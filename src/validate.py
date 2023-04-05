# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 12:02:27 2023

@author: ZimmermannP
"""

"""verschiedene experimentiell bestimmte Wärmeübergänge von Ammoniak
   microfin/plain tube aus Diagrammen abgelesen von verschiedenen Autoren"""
import matplotlib.pyplot as plt
from src.model_zürcher import alpha_m

def zürcher_m10():
    """ Zürcher Diss plain tube
    di = 0.014
    s = 1/1000
    te = 4
    q = 7130
    m = 10
    lw = 26"""
   
    x = [0.28,0.31,0.35,0.47,0.5,0.6,0.8]
    a = [1400,1350,1270,1200,1170,1100,1000]
    
    if plot == True:
        val =  alpha_m(di, s, lw, te, m, q, t)
        plt.plot(val[:,0], val[:,1], label = (f'Korrelation Zürcher {m} [kg/m²s] {round(q/1000,2)} [kW/m²] Glattrohr'))
        
    return x,a
    
def zürcher_m20_5400():
    """ Zürcher Diss plain tube
    di = 0.014
    te = 4
    q = 5400
    m = 20
    lw = 26"""
    x = [0.07,0.11,0.13,0.145,0.175,0.176,0.19,0.23,0.32,0.345,0.36,0.39,0.41,0.43,0.45,0.51,0.52,0.59,0.68,0.73,0.82,0.85,0.95]
    a = [1700,1900,1800,1150,1200,1850,1300,1450,940,1000,700,900,700,650,750,720,680,700,600,720,515,490,500]
    return x,a
def zürcher_m20_1780():
    """ Zürcher Diss plain tube
    di = 0.014
    te = 4
    q = 17800
    m = 20
    lw = 26"""
    x = [0.345,0.36,0.395,0.45,0.48,0.495,0.51,0.52,0.62,0.63,0.71,0.725,0.775,0.795,0.805,0.86,0.88,0.890,0.96]
    a = [2990,2920,2780,2750,2760,2700,2600,2550,2500,2430,2300,2360,2100,1900,1950,1750,1550,1650,1500]
    return x,a

def validate_zürcher():
    """ plot zur validierung von Zürcher """
    di = 0.014
    s = 1/1000
    te = 4
    lw = 26
    q1 = 17800
    q2 = 5400
    q3 = 7130
    m1 = 20
    m2 = 10
    t = 50

    val =  alpha_m(di, s, lw, te, m1, q1, t)
    val2 = alpha_m(di, s, lw, te, m1, q2, t)
    val3 = alpha_m(di, s, lw, te, m2, q3, t)

    xd1,alpha_z1 = zürcher_m10()
    xd2, alpha_z2 = zürcher_m20_5400()
    xd3,alpha_z3 = zürcher_m20_1780()

    plt.plot(val[:,0], val[:,1], label = (f'Korrelation Zürcher {m1} [kg/m²s] {round(q1/1000,2)} [kW/m²] Glattrohr'))
    plt.plot(val[:,0], val2[:,1], label = (f'Korrelation Zürcher {m1} [kg/m²s] {round(q2/1000,2)} [kW/m²] Glattrohr'))
    plt.plot(val[:,0], val3[:,1], label = (f'Korrelation Zürcher {m2} [kg/m²s] {round(q3/1000,2)} [kW/m²] Glattrohr'))
    plt.scatter(xd3,alpha_z3,label = (f'Messwerte Zürcher {m1} [kg/m²s] {round(q1/1000,2)} [kW/m²] Glattrohr'))
    plt.scatter(xd2,alpha_z2,label = (f'Messwerte Zürcher {m1} [kg/m²s] {round(q2/1000,2)} [kW/m²] Glattrohr'))
    plt.scatter(xd1,alpha_z1,label = (f'Messwerte Zürcher {m2} [kg/m²s] {round(q3/1000,2)} [kW/m²] Glattrohr'))


    plt.ylim(0)
    plt.xlim(0,1)
    plt.ylabel('\u03B1_i [W/m²K]')
    plt.xlabel('Massendampfgehalt [-]')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
    #plt.title(f'{R}; {te}°C;')
    plt.title(f'R717, Tsat={te}[°C], di={round(di*1000)}[mm] ')
    plt.grid(True)
    
def kabelac_m50(plain=True):
    """ te = -20
        m = 50
        q = 40 kW
        di = 11.3 mm
        drallwinkel y = 25
        25 Rippen
        höhe Rippen = 0.63mm """
    if plain == False:
        x = [0.16,0.35,0.52,0.9]
        a = [12000,14500,16000,13000]
    else: 
        x = [0.08,0.12,0.18,0.21,0.28,0.325,0.36,0.42,0.46,0.5,0.58,0.61,0.69,0.78,0.79,0.82,0.88,0.9,0.96]
        a = [4100,4500,4900,5000,6000,6200,7000,7400,8000,8800,9900,10000,10200,10500,11000,11200,9000,8000,3500]
    return x,a
def kelly_m27_20():
    """ Werte nach Kelly,Eckels,Fenton 2002
    Alurohr
        te = -20
        m = 27
        q = 2710 """
    x = [0.3,0.475,0.74,0.755]
    a = [3250,5600,7250,7000]
    return x,a

def kelly_m27_10():
    x = [0.12,0.225,0.51,0.78,0.87]
    a = [2800,3800,6750,5700,5100]
    return x,a
    
def kelly_m9_10():
    """ Werte nach Kelly,Eckels,Fenton 2002 """
    """ te = -10
        m = 9
        q = 2710 """
    x = [0.41,0.675]
    a = [3500,3750]
    return x,a
    
# def validate_kabelac():
        
        
        


