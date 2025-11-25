import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

# %%
# Set the font to Times and enable LaTeX rendering with pdflatex
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times']
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{newtxmath}'
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['pgf.texsystem'] = 'pdflatex'

plt.rcParams["axes.prop_cycle"] = cycler(color="brgcmyk")

# %%
# Calculate the electrical parameters (constants)

T = 45  # Temperature in Celsius degrees
N = 22  # Number of cells in series

# Ohmic resistance
A = 300e-4  # Cell surface in m^2
r1 = 59.5482e-6
r2 = -340.8224e-9
r3 = -106.9708e-6
r4 = 2.7075e-3
R_ohm = N*(r1 + r2*T + r3/T + r4/T**2)/A

# Anode
C_an = .6375
s_an1 = 25.23e-3 
#s_an1 = 1.64*25.23e-3  # Correction corresponding to Fig. 11 in Ursua and Sanchis (2011)
s_an2 = -234.0338e-6
s_an3 = 3.1832e-6
t_an1 = 54.6185e-3
t_an2 = -2.4601e-3
t_an3 = 52.1217e-6
# Calculate the incremental activation resistance at zero current, R_an0 = R_an(0)
J_an = t_an1 + t_an2*T + t_an3*T**2 # Current at which R_an = 0.5*R_an0
V_an = N*(s_an1 + s_an2*T + s_an3*T**2)
R_an0 = V_an/J_an

# Cathode activation
C_ca = .0307
s_ca1 = 110.3623e-3 
#s_ca1 = 1.13*110.3623e-3  # Correction for Fig. 11
s_ca2 = -1.6466e-3
s_ca3 = 22.8382e-6
t_ca1 = 45.7027 
#t_ca1 = 1.2028*45.7027  # Correction for Fig. 11
t_ca2 = .7781
t_ca3 = -10.5743e-3
# Calculate the incremental activation resistance at zero current, R_ca0 = R_ca(0)
J_ca = t_ca1 + t_ca2*T + t_ca3*T**2  # Current at which R_ca = 0.5*R_ca0
V_ca = N*(s_ca1 + s_ca2*T + s_ca3*T**2)
R_ca0 = V_ca/J_ca


# %%

def R_an(i):
    """Incremental anode activation resistance at the operating point."""
    return R_an0/(1 + i/J_an)

def R_ca(i):
    """Incremental cathode activation resistance at the operating point."""
    return R_ca0/(1 + i/J_ca)

i_list = [20, 70, 120]

w = 2*np.pi*np.logspace(-4, 5, 1000)
s = 1j*w

Z = []
for i in i_list:
    Z.append(R_ohm + R_an(i)/(s*R_an(i)*C_an + 1) + R_ca(i)/(s*R_ca(i)*C_ca + 1))


# %%
#plt.figure(figsize=(3, 2.5))
for index, i in enumerate(i_list):
    plt.plot(Z[index].real, Z[index].imag, label=f'{i} A', lw=1.5)

plt.xlabel(r'$\mathrm{Re}\{ Z(\mathrm{j}\omega) \} \ (\Omega)$')
plt.ylabel(r'$\mathrm{Im}\{ Z(\mathrm{j}\omega) \} \ (\Omega)$')
# plt.legend(loc='center left')
plt.legend(loc='best')
plt.grid()
plt.axis('image')
# plt.axis([0, .1, -.03, 0])
# plt.yticks([0, -.02])
plt.axis([0, .08, -.02, 0])
plt.yticks([0, -.01, -.02])
plt.tight_layout()
plt.show()

#plt.savefig('nyqs.pdf')

# %%
# For checking purposes
# Values from Ursua and Sanchis (2011), Fig. 11
i_list = [20, 70 ,120]
R_an_list = [.041, .0135, .0085]
R_ca_list = [.024, .0145, .0127]

i_values = np.linspace(20, 100, 500)
# plt.figure(figsize=(8, 6))
plt.plot(i_values, R_ca(i_values), label='R_ca(i)', lw=2)
plt.plot(i_list, R_ca_list, label='R_ca_list', marker='o', linestyle='')
plt.axis([0, 120, 0, .05])
plt.xlabel('Current (A)')
plt.ylabel('R_ca (Ohm)')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

plt.plot(i_values, R_an(i_values), label='R_an(i)', lw=1.5)
plt.plot(i_list, R_an_list, label='R_an_list', marker='o', linestyle='')  
plt.axis([0, 120, 0, .05])
plt.xlabel('Current (A)')
plt.ylabel('R_an (Ohm)')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()