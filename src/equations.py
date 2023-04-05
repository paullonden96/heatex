# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 14:58:28 2022

Formeln aus Plan
"""

    
    
    
    
    '''
    #latente Wärmestrom am Grundrohr bezogen auf die Reifoberfläche A_F_G
    def Q_l_G(t_G_m):
        return gr.F3_G * air.beta * ((air.p_D_L() - gr.p_D_F_G_quer)\
                / (air.R_D*(273.15+(air.v_L+gr.t_f_G_quer/2)))*air.delta_h_sub)
            
    #sensible Wärmestrom am Grundrohr bezogen auf die Reifoberfläche
    def Q_s_G():
        F4 = 2.5*(air.lambda_L/(geo.alpha_a()*gr.sigma_f_g))
        alpha_f = F4*geo.alpha_a()
        return alpha_f * (air.v_L - (air.t_F_o(gr.t_G) + gr.t_G)/2)
        
   '''
"""        
'Rippen'   
class r:
    '''
    #t_R = -20.5 #mittlere Rippentemperatur Schätzung
    t_R = air.v_L-(air.v_L-gr.t_G)*geo.eta_R(alpha_a) #mittlere Rippentemperatur über Rippenwirkungsgrad aktuell -24°C
    sigma_f_r = air.sigma_f(t,t_R) #Reifdicke an den Rippen
    t_F_o_R = air.t_F_o(t_R) #Temperatur an der Reifspitze für Rippen
    p_D_R_strich = air.p_D_strich(t_R) #Sättigungsdruck bei t_R an den Rippen
    p_D_F_R_quer = air.p_D_F_quer(t_R) #wirksame Wasserdampfdruck im Reif Rippen
    t_F_R_quer = air.t_f(p_D_F_R_quer) #Sättigungstemp. bei p_D_F_R_quer
    F3_R = 1.8*(air.D / (air.beta*sigma_f_r))**0.5 #Faktor F3 Rippen
    '''
    
    #latente Wärmestrom am Grundrohr bezogen auf die Reifoberfläche A_F_R
    def Q_l_R():
        return r.F3_R * air.beta * ((air.p_D_L() - r.p_D_F_R_quer)\
                / (air.R_D*(273.15+(air.v_L+gr.t_f_G_quer/2)))*air.delta_h_sub)
            
    #sensible Wärmestrom am Grundrohr bezogen auf die Reifoberfläche A_F_R
    def Q_s_R(alpha_a):
        F4 = 2.5*(air.lambda_L/(alpha_a*r.sigma_f_r))
        alpha_f = F4*alpha_a
        return alpha_f * (air.v_L - (air.t_F_o(r.t_R) + r.t_R)/2)

    

 
def tlm(t1,tr,k,A,M,dhdl,p):
    t_l_m = t1 - (t1-tr) * (1-math.exp(-((k*A)/(p*M*dhdl))))
    return t_l_m

if __name__ == "__main__":
         
"""