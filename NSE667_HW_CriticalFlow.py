# script for NSE 667 Critical Flow HW
# Austin Warren
# June 2022

import numpy as np
import matplotlib.pyplot as plt
from pyXSteam.XSteam import XSteam
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)

def MassFlux(h_0,h_m,rho_g,rho_f,x,S):
    rho_m = (((x/rho_g)+((1-x)/rho_f)*S)*(x + (1-x)*(1/S**2))**(1/2))**(-1)
    G_m = rho_m * np.sqrt(2*(h_0 - h_m)*1000)
    return G_m



# inital conditions
P_0 = 35; # initial pressure in bar
T_0 = steamTable.tsat_p(35) # initial temperature in deg C
h_0 = steamTable.hL_p(35) # initial saturated enthlapy in kj/kg
s_0 = steamTable.sL_p(35) # initial specific entropy


# generate pressure values from 0 to P_0
# use pressure values and entropy to determine corresponding h_m
# calc G_m for the pressure
P = np.linspace(0.01,P_0, 100)
x = np.zeros(len(P))

# calc rho_m
rho_g = np.zeros(len(P))
rho_f = np.zeros(len(P))
rho_m = np.zeros(len(P))

# calc h_m
h_g = np.zeros(len(P))
h_f = np.zeros(len(P))
h_m = np.zeros(len(P))
for i in range(len(P)):
    x[i] = steamTable.x_ps(P[i],s_0)

    rho_g[i] = steamTable.rhoV_p(P[i])
    rho_f[i] = steamTable.rhoL_p(P[i])
    rho_m[i] = (x[i]/rho_g[i] + (1-x[i])/rho_f[i])**(-1)

    h_g[i] = steamTable.hV_p(P[i])
    h_f[i] = steamTable.hL_p(P[i])
    h_m[i] = x[i]*h_g[i] + (1-x[i])*h_f[i]

# G_m -- HEM
G_m = rho_m * np.sqrt(2*(h_0 - h_m)*1000)


# Fauske
S_fauske = np.sqrt(rho_f / rho_g)
G_m_fauske = MassFlux(h_0,h_m,rho_g,rho_f,x,S_fauske)

# Moody
S_moody = (rho_f / rho_g)**(1/3)
G_m_moody = MassFlux(h_0,h_m,rho_g,rho_f,x,S_moody)



# plot
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.figure(facecolor='w', edgecolor='k', dpi=300)
plt.plot(G_m, P/P_0, '-k', label='HEM')
plt.plot(G_m_fauske, P/P_0, '-b', label='Fauske')
plt.plot(G_m_moody, P/P_0, '-r', label='Moody')

plt.xlabel(r'$G_m$ ($\frac{kg}{m^2s}$)')
plt.ylabel(r'$\frac{P}{P_0}$')
plt.figlegend(bbox_to_anchor=(1.0,0.9))
plt.grid(b=True, which='major', axis='both')
plt.savefig('1-plots/graph_Gm.pdf',transparent=True)