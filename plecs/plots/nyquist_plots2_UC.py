import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
import sympy as sp
from scipy import signal

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
T = 45
N_s = 22
A = 300e-4

s1a = 25.23e-3
s2a = -234.0338e-6
s3a = 3.1832e-6
t1a = 54.6185e-3
t2a = -2.4601e-3
t3a = 52.1217e-6

s1c = 110.3623e-3
s2c = -1.6466e-3
s3c = 22.8382e-6
t1c = 45.7027
t2c = .7781
t3c = -10.5743e-3

# %%
def parameter_s(T, s1, s2, s3):
    return s1 + s2*T + s3**T**2

def parameter_t(T, t1, t2, t3):
     return t1 + t2*T + t3*T**2

def activation_overvoltage(T, s1, s2, s3, t1, t2, t3, i):
    return parameter_s(T, s1, s2, s3)*np.log((i/parameter_t(T, t1, t2, t3)) + 1)

def activation_current(T, s1, s2, s3, t1, t2, t3, u_act):
    return parameter_t(T, t1, t2, t3)*(np.exp(u_act/(parameter_s(T, s1, s2, s3))) - 1)
##############
# This gives the so called the chord-slope resistance, i.e., the straight line that connects the origin to the operating point. It should be the incremental resistance at the operating point instead (the slope). 
# def activation_resistance(T, s1, s2, s3, t1, t2, t3, i, N_s):
#     u_act = activation_overvoltage(T, s1, s2, s3, t1, t2, t3, i, N_s)
#     i_act = activation_current(T, s1, s2, s3, t1, t2, t3, i, N_s, u_act)
#     print(i_act)
#     return u_act/i_act

def activation_resistance(T, s1, s2, s3, t1, t2, t3, i, N_s):
    u_act = activation_overvoltage(T, s1, s2, s3, t1, t2, t3, i)
    i_act = activation_current(T, s1, s2, s3, t1, t2, t3, u_act)
    return N_s * (parameter_s(T, s1, s2, s3)/ (i_act + parameter_t(T, t1, t2, t3)))

def ohmic_resistance(T, N_s, A):
    r1 = 59.5482e-6
    r2 = -340.8224e-9
    r3 = -106.9708e-6
    r4 = 2.7075e-3
    r = r1 + r2*T + r3/T + r4/T**2
    return N_s*r/A

#i_list = [40, 80, 120]
i_list = [120] #, 120]

#R_ohm = ohmic_resistance(T, N_s, A)
R_ohm = 0.032
C_an = .6375
C_ca = .0307

w = 2*np.pi*np.logspace(-4, 5, 1000)
s = 1j*w
Z = []

for i in i_list:
    #R_an = activation_resistance(T, s1a, s2a, s3a, t1a, t2a, t3a, i, N_s)
    #R_ca = activation_resistance(T, s1c, s2c, s3c, t1c, t2c, t3c, i, N_s)
    R_an = 0.023
    R_ca = 0.023
    Z.append(R_ohm + R_an/(s*R_an*C_an + 1) + R_ca/(s*R_ca*C_ca + 1))



# %%
plt.figure(figsize=(3, 2.5))

for index, i in enumerate(i_list):
   plt.plot(w, abs(Z[index]), label=f'{i} A', lw=2)
#plt.semilogx(w_g, mag)

plt.xlabel(r'$\mathrm{Re}\{ Z(\mathrm{j}\omega) \} \ (\Omega)$')
plt.ylabel(r'$\mathrm{Im}\{ Z(\mathrm{j}\omega) \} \ (\Omega)$')
#plt.legend(loc='lower left')
plt.grid()
plt.axis('image')
#plt.axis([0, .16, -.08, 0])
#plt.xticks([0, .01, .02, .03, .04, .05, .06, .07, .08])
#plt.yticks([0, -.004, -.008, -.012, -.016])
# plt.xlim([0.03, 0.14])
#plt.xticks(np.arange(0.03, 0.14, 0.01))
# plt.ylim([-0.040, 0.000])
plt.tight_layout()
plt.show()
#plt.savefig('nyqs.pdf')
#print(R_ohm)

