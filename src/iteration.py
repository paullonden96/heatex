# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 09:23:29 2022

TO-DO S.489 das Beispiel durchrechnen mit t=0 beginnen und t=1,2,3 herausfinden 
Ziel Tabelle wie im BSP danach inneres alpha!

Ziel der Arbeit:
- Berechnungen zu Strömungsgeschwindigkeiten, Druckverlusten, Wärmeübergängen 


@author: ZimmermannP
"""

from src.air import x, air

from src.model_zürcher import alpha_m
from src.deltap import deltap_b, deltap_r_G
from src.geometrie import etaR
from src.equations import t_R_m,tg_m,alpha_scheinbar,k,tm



def it_tl(tl,k_s,alpha_as,t_r_m,alpha_a,Q_h,eta_R,air1,geo1,alpha_i,t_strich_l,Aa_s):
    """ Iteration der mittleren Lufttemperatur im Verdampfer nach Plank Kap 12. S. 482
    Schritt 13 - 21 """

    t_g_m = tg_m(tl,k_s,alpha_as,t_r_m)
    print(f'neue Grundrohrtempertur {t_g_m}')
    #13-16
    alpha_G = air1.alpha_g_s(Q_h,t_g_m,alpha_a)
    #17
    t_r_m = t_R_m(tl, t_g_m, eta_R)
    print(f'neue Rippentemp. {t_r_m}')
    #18-20
    alpha_R = air1.alpha_g_s(Q_h,t_r_m,alpha_a)
    #21
    eta_R = etaR(geo1,alpha_R)
    
    dhdt_s = air1.dhdt(Q_h,geo1,t_r_m,t_g_m)
    print(dhdt_s)
    alpha_as = alpha_scheinbar(alpha_G,geo1,alpha_R,eta_R)
    print('alpha_as=',alpha_as)
    k_s = k(alpha_as,alpha_i,geo1)
    
    t_l_m = tm(t_strich_l,t_r_m,k_s,Aa_s/2,air1.M_L,dhdt_s)
    #print(f'mittlere Lufttemp. {i}.Durchlauf {t_l_m}')
    """ 2.Durchlauf """
    t_g_m = tg_m(t_l_m,k_s,alpha_as,t_r_m)
    print(f'neue Grundrohrtempertur {t_g_m}')
    
    alpha_G = air1.alpha_g_s(Q_h,t_g_m,alpha_a)
    #17
    t_r_m = t_R_m(t_l_m, t_g_m, eta_R)
    print(f'neue Rippentemp. {t_r_m}')
    #18-20
    alpha_R = air1.alpha_g_s(Q_h,t_r_m,alpha_a)
    #21
    eta_R = etaR(geo1,alpha_R)
    
    dhdt_s = air1.dhdt(Q_h,geo1,t_r_m,t_g_m)
    
    alpha_as = alpha_scheinbar(alpha_G,geo1,alpha_R,eta_R)
    print('alpha_as=',alpha_as)
    k_s = k(alpha_as,alpha_i,geo1)
    
    t_l_m = tm(t_strich_l,t_r_m,k_s,Aa_s/2,air1.M_L,dhdt_s)
    print(f'mittlere Lufttemp. 2.Durchlauf {t_l_m}')
    
    
        
    
    



