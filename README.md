# heatex
flowboiling thermodynamic calculations for an ventilation evaporator:
  - flow pattern map by Steiner (VDI) and Zürcher (dissertation [1])
  - heat transfer coefficient flow boiling calculated by Steiner (VDI Wärmeatlas) and Zürcher [1]
  - heat transfer for air (forced convektion) by Plank [2]
  - heat transfer flow boiling microfinned by Chamra [3] and Thome [4]

# calculation of alpha_i (heat transfer inside) 
``` python
""" Example: """

di = 0.012  # inside diameter [m]
s = 1 / 10**3  # thickness of pipe [m]
m = 15  # mass velocity [kg/m²s]
lw = 16 # conductivity pipe [W/m²K]
te = -30 # evaporation temperature [°C]
q = 6400 # heat flux [W/m²]
t = 50 # steps

# val is a matrix you get from calling alpha_m (Zürcher) or alpha_m_s (Steiner)
# alpha_m calculates for every x (vapor quality) the heat transfer inside and gives you a matrix with each pair x-alpha

val = alpha_m(di, s, lw, te, m, q, t)

# the plot
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
```
# calculation of heat transfer outside 
## making a geometrie object of the evaporator 
```python
t_R = 0.004 # Lamellenteilung [m]
sigma_R = 0.00025 # Lamellendicke [m]
di = 0.0155 # Innendurchmesser Grundrohr [m]
da = 0.0165 # Außendurchmesser Grundrohr [m]
s = 0.001 # Wandstärke Grundrohr [m]
sq = 0.05 # Abstand Rohrachsen Verdampferrohre senkrecht [m]
sl = 0.05 # Abstand Rohrachsen Verdampferrohre waagerecht [m]
i_R = 12 # Anzahl Rohrreihen
i_Q = 4 # Anzahl Rohre pro Rohrreihe
z = 6 # parallelgeschaltete Rohreihen
Al = 1.2 # Fläche Ventilatorkanal [m²]
R_p = 10*-6 # technisch glatte Rohre?
L = 2 # Rohrlänge [m] --> Teilsegmente für die alpha berechnet werden soll
lR = 220 # Wärmeleitfähigkeit Rippe --> Alu [W/m²K]
lG = 16 # Wärmeleitfähigkeit Grundrohr --> Stahl [W/m²K]
m = 15 # Massenstromdichte [kg/m²s]
tg = te # Grundrohrtemperatur tg = te vereinfacht angenommen
tl = -18 # air temperatur entry [°C]
phi_i = 0.91 # relative air humidity [-]
Vl = 10000 # air flow rate [m³/h]


