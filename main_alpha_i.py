# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 11:42:20 2022

@author: ZimmermannP

Berechnungsprogramm für den inneren Wärmeübergang nach verschiedenen Modellen 
und unter Einfluss verschiedener Rohroptimierungen 

Voraussetzung: 
    pip install coolprop
    pip install numpy
    pip install matplotlib

Bei Fehlermeldung "No Module named 'src' " in Spyder muss der Ordner 'heatex' geöffnet werden.
Ganz oben rechts auf das Symbol was wie eine Akte aussieht klicken und den Ordner auswählen""
Empfehlung: 
    Nutzung von Spyder zum automatischen Anzeigen der plots
    https://www.spyder-ide.org/
    
"""

from src.Refrigerant import refrigerant 
from src.model_zürcher import strömungskarte, alpha_m
from src.model_steiner import strömungskarte_s, alpha_m_s
from src.deltap import deltap_r_m
from src.model_chamra_microfin import alpah_tf_m, Kabelac, Kelly
from src.model_thome_microfin import alpha_tf_m_thome
from src.validate import *


import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import math
import numpy as np
Pi = math.pi



''' Berechnung eines mittleren alpha_i bei gleicher Massenstromdichte
    alpha, val = alpha_m(di,s,lw,*params,ausgebildet)
    Rückgabe von alpha gemittelt und val Matrix Dimension (t,11)
    val[:,0] = x Dampfmassengehalt
    val[:,1] = Stelle z im Rohr 
    val[:,2] = alpha_z lokaler Wärmeübergangskoeffizient an der Stelle z ACHTUNG: Annahme linearer Anstieg von x über z
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

# L = 3
# di = 0.014
# s = 0.93/1000
# lw = 28
# te = 4
# m = 20


# xmax = 1 # x nach Verdampfung
# t = 100 #Anzahl der Abschnitte in die das Rohr unterteilt werden soll
# lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l = ref.stoffwerte(te, R)

# ausgebildet = False # Berechnung mit nicht ausgebildeter Strömung siehe VDI H3.5 Anhang
# params = [t,L,R,Ra,theta,q,lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe]
# alpha1, val1 = alpha_m_s(di,s,lw,te,m,*params)
# xd = val1[:,0]
# dpdl = val1[:,9]
# alpha_z = val1[:,2]
# alpha_z_k_l0 = val1[:,5]

# print(f'\u03B1_i_m = {alpha1} W/m²K --> Glattrohr')
# print(f'\u03B1_i_m = {alpha2} W/m²K --> Twisted Tape')


''' print alpha(z)_k / alpha_L0 Verhältnis von alpha(z)/alpha_l0 '''
# ylabel = '\u03B1(z)_k/\u03B1_L0'
# plt.plot(xd,alpha_z_k_l0)
# plt.ylabel('\u03B1(z)_k/\u03B1_L0')
# plt.xlabel('Massendampfgehalt x')
# plt.ylim(0.5,20)


'''print alpha(Z)'''
# if val1[0,3] != 0:  plt.plot(xd,val1[:,3], label = '\u03B1(Z)_B')
# plt.plot(xd,alpha_z)
# #plt.axhline(y=alpha1, color = 'red', label = '\u03B1_i_m' )
# #plt.axhline(y=alpha2, color = 'green', label = '\u03B1_i_m')
# plt.ylabel('\u03B1 in W/m²K')
# plt.xlabel('Massendampfgehalt x')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
# plt.grid(True)
# plt.title(f' {R}; {m}kg/m²s; {te}°C; {round(di*1000,2)}mm; {L}m')


''' Berechnung mittleren alphas für alle Massenstromdichten von m1 - m2 '''
# ausgebildet = False # Berechnung mit nicht ausgebildeter Strömung siehe VDI H3.5 Anhang
# m1 = 1 
# m2 = 30
# mx_16_al,ax_16_al = alpham_m(m1,m2,0.016,0.0005,16,*params,ausgebildet) # 16mm Innndurchmesser Alulegierung 220 W/m²K
# # mx_12_al,ax_12_al = alpham_m(m1,m2,0.012,0.0015,220,*params,ausgebildet)

# # mx_16_st,ax_16_st = alpham_m(m1,m2,0.016,0.0005,16,*params,ausgebildet) # 16mm Innendurchmesser Stahlrohr 
# # mx_12_st,ax_12_st = alpham_m(m1,m2,0.012,0.0005,16,*params,ausgebildet)

# plt.plot(mx_16_al,ax_16_al, label=(f'Alu {round(0.016*1000,2)}x{0.0015*1000}mm; {220}W/m²K; n.ausg.'))
# # plt.plot(mx_12_al,ax_12_al, label=(f'Alu {round(0.012*1000,2)}x{0.0012*1000}mm; {220}W/m²K; n.ausg.'))

# # plt.plot(mx_12_st,ax_12_st, label=(f'Stahl {round(0.012*1000,2)}x{0.0005*1000}mm; {16}W/m²K; n.ausg.'))
# # plt.plot(mx_16_st,ax_16_st, label=(f'Stahl {round(0.016*1000,2)}x{0.0005*1000}mm; {16}W/m²K; n.ausg.'))
# plt.title(f'{R}; {te}°C; {L}m Rohrlänge')
# plt.ylabel('\u03B1_i_m in W/m²K')
# plt.xlabel('Massenstromdichte in kg/m²s')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
# plt.grid(True)


''' Druckverlust in Abhängigkeit von x '''
# m = 15
# te = -30
# t = 30
# di1 = 0.008 # Meter
# di2 = 0.010
# di3 = 0.012
# di4 = 0.014

# val1 = deltap_r_m(m,te,di1,t) #Glattrohr, 16mm ausgebildet = False
# val2 = deltap_r_m(m,te,di2,t) 
# val3 = deltap_r_m(m,te,di3,t) 
# val4 = deltap_r_m(m,te,di4,t) #Glattrohr, 12mm ausgebildet = False
# xd = val1[:,0]
# dpdl1 = val1[:,1] / 100
# dpdl2 = val2[:,1] / 100
# dpdl3 = val3[:,1] / 100
# dpdl4 = val4[:,1] / 100

# plt.plot(xd,dpdl1, label = (f'Stahl {di1*1000:.2f} x {0.0005*1000}mm Glattrohr'))
# plt.plot(xd,dpdl2, label = (f'Stahl {di2*1000:.2f} x {0.0005*1000}mm Glattrohr'))
# plt.plot(xd,dpdl3, label = (f'Stahl {di3*1000:.2f} x {0.0005*1000}mm Glattrohr'))
# plt.plot(xd,dpdl4, label = (f'Stahl {di4*1000:.2f} x {0.0005*1000}mm Glattrohr'))
# plt.ylabel('\u0394p / \u0394L [mbar/m]')
# plt.xlabel('Massendampfgehalt [-]')
# plt.title(f'R717; {te}°C; {m}kg/m²s')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
# plt.grid(True)
# plt.xlim(0,1)
# plt.ylim(0)


''' Vergleich von benetzt und unbenetzt '''
# di1 = 0.0155
# s1 = 0.0005
# m1 = 15
# m2 = 40
# te = -30

# alpham1, val1 = alpha_m(di1,s1,28,te,m1,*params) 
# xd = val1[:,0]
# zd = val1[:,1]
# az = val1[:,2]


# alpham2, val2 = alpha_m(di1,s1,28,te,m2,*params)
# xd1 = val2[:,0]
# zd1 = val2[:,1]
# az1 = val2[:,2]

# plt.plot(xd,az, label = (f'Stahl {m1} [kg/m²s] {di1*1000:.2f} x {s1*1000}mm Glattrohr'))
# plt.plot(xd,az1, label =(f'Stahl {m2} [kg/m²s] {di1*1000:.2f} x {s1*1000}mm Glattrohr'))
# plt.ylabel('\u03B1_i [W/m²K]')
# plt.xlabel('Massendampfgehalt [-]')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
# plt.title(f'{R}; {te}°C;')
# plt.grid(True)

''' Strömungskarte nach Zürcher '''

# di = 0.0155
# te = -30
# m = 15
# q = 6500
# t = 50 

# lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv = ref.stoffwerte(te, R)
# params = [t,R,Ra,theta,q,lambda_g,lambda_l,rho_l,rho_g,eta_g,eta_l,pr_g,sigma_r,cp_g,cp_l,delta_hv,pe]
# val, x = strömungskarte(di,0.0005,16,te,m,*params)
# xd = val[:,0]
# zd = val[:,1]
# G_strat = val[:,2]
# G_wavy = val[:,3]

# G_wavy_min = min(val[:,3]) / 180
# ring_text_y = (1-G_wavy_min) - G_wavy_min
# fontsize = 'large'
# #plt.plot(xd,m_s, label = 'Grenzkurve Schichtenströmung')
# plt.plot(xd,G_strat,label = 'Grenzkurve Schichten/Wellenströmung Zürcher')
# plt.plot(xd,G_wavy, label = 'Grenzkurve Wellen/Ringströmung Zürcher')
# #plt.plot(xd,G_wavy_ohne_dis, label = 'Grenzkurve Wellen/Ringströmung Zürcher ohne Dissipation')
# plt.figtext(0.15,0.135,'Schichtenströmung',fontsize=fontsize)
# plt.figtext(0.35,G_wavy_min,'Schichten-Wellenströmung',fontsize=fontsize)
# plt.figtext(0.4,ring_text_y,'Ringströmung',fontsize=fontsize)

# #plt.plot(xd,m_w_s, label = 'Grenzkurve Wellenströmung nach Steiner', color = 'black')
# #plt.plot(xd,m_w_k, label = 'Grenzkurve Wellenströmung nach Kattan')
# #plt.plot(xd,m_w_z, 'r--', label = 'Grenzkurve Wellenströmung nach Zürcher')
# plt.vlines(x= x[0], ymin=x[1], ymax=180, label ='Grenzkurve Ringströmung', color = 'green') 
# plt.ylabel('Massenstromdichte [kg/m²s]')
# plt.xlabel('Massendampfgehalt [-]')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
# plt.title(f'{R}; te={te}°C; di={di*1000}mm, q={q/1000} kW/m² Glattrohr')
# plt.ylim(0,100)
# plt.xlim(0,1)
# plt.grid(True)

""" Strömungskarte nach Zürcher/Steiner """
# di = 0.0118
# s = 1/1000
# te = -20
# m = 9
# q = 2710    
# t = 50
# lw = 220
# L = 1
# val,phi_wet, x_A, G_wavy_A = strömungskarte_s(di,s,lw,te,m,q,t,L)
# val2,x = strömungskarte(di,s,lw,te,m,q,t)
# xd = val[:,0]
# G_strat = val[:,1]
# G_wavy = val[:,2]
# xd = val2[:,0]
# G_wavy_Z = val2[:,3]
# a_T = val2[:,4]
# a_R = val2[:,5]
# G_wavy_min = min(val[:,2]) / 180
# ring_text_y = (1-G_wavy_min) - G_wavy_min
# fontsize = 'large'

# #plt.plot(xd,a_R,label = 'Dampfvolumen Rouhani')
# #plt.plot(xd,m_s, label = 'Grenzkurve Schichtenströmung')
# plt.plot(xd,G_strat,label = 'Grenzkurve Schichten/Wellenströmung')
# plt.plot(xd,G_wavy, label = 'Grenzkurve Wellen/Ringströmung Steiner VDI')
# plt.plot(xd,G_wavy_Z,'--', color='r', label = 'Grenzkurve Wellen/Ringströmung Zürcher')

# plt.figtext(0.15,0.135,'Schichtenströmung',fontsize=fontsize)
# plt.figtext(0.35,G_wavy_min,'Schichten-Wellenströmung',fontsize=fontsize)
# plt.figtext(0.4,ring_text_y,'Ringströmung',fontsize=fontsize)

# plt.vlines(x= x_A, ymin=G_wavy_A, ymax=200, label ='Grenzkurve Ringströmung', color = 'green') 
# plt.ylabel('Massenstromdichte [kg/m²s]')
# plt.xlabel('Massendampfgehalt [-]')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
# plt.title(f'R717, te={te}°C; di={di*1000}mm, q={q/1000}kW/m² Glattrohr')
# plt.ylim(0,160)
# plt.xlim(0,1)
# plt.grid(True)

''' alpha nach Zürcher '''
di = 0.012
s = 1 / 10**3
m = 15
lw = 16
te = -30
q = 6400
t = 50
val = alpha_m(di, s, lw, te, m, q, t)

xd = val[:,0]
alpha_z = val[:,1]
a_m = np.mean(val[:,1])
print(a_m)
plt.plot(val[:,0], val[:,1], label = (f'Stahl {lw} [W/mK], {m} [kg/m²s] {di*1000:.2f} x {s*1000} [mm] Glattrohr'))
plt.ylim(0,max(alpha_z)+100)
plt.xlim(0,1)
plt.ylabel('\u03B1_i [W/m²K]')
plt.xlabel('Massendampfgehalt [-]')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
plt.title(f'R717; t_sat={te}°C; q={q/1000}kW/m² Glattrohr')
plt.grid(True)

""" Vergleich  Rippenrohr Experiment/Theorie Messwerte Kabelac """
# te = -20
# geo = Kabelac()
# di = 0.01113
# s = 1/1000
# q = 40000
# m = 50
# t = 20
# val = alpah_tf_m(m,q,di,s,te,t,Kabelac())
# val2 = alpha_tf_m_thome(m,q,di,s,te,t,Kabelac())
# xd = val[:,0]
# alpha_z = val[:,1]
# xd2, alpha_z2 = kabelac_m50(plain=False)
# xd3,alpha_z3 = kabelac_m50(plain=True)
# alpha_z4 = val2[:,1]

# plt.plot(val[:,0], val[:,1], label = ('berechnet Rippenrohr nach Chamra'))
# plt.plot(val[:,0], val2[:,1], label = ('berechnet Rippenrohr nach Thome'))
# plt.scatter(xd2,alpha_z2,label = ('Messwerte Rippenrohr Kabelac'))
# plt.scatter(xd3,alpha_z3,label = ('Messwerte Glattrohr Kabelac'))
# plt.ylim(0,max(val2[:,1])*1.1)
# plt.xlim(0,1)
# plt.ylabel('\u03B1_i [W/m²K]')
# plt.xlabel('Massendampfgehalt [-]')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) 
# plt.title(f'R717,Tsat= {te}[°C], {m}[kg/m²s], {q/1000} [kW/m²], Vergleich Korrelation/Messwerte Kabelac')
# plt.grid(True)

""" Validierung Modell Zürcher mit seinen Messwerten """

# validate_zürcher() # plot von Graphen zur Validierung von Zürcher

