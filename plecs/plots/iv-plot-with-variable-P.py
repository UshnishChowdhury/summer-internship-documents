import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
import math

# %%
# Set the font to Times and enable LaTeX rendering with pdflatex
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times']
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{newtxmath}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{gensymb}'
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['pgf.texsystem'] = 'pdflatex'

plt.rcParams["axes.prop_cycle"] = cycler(color="brgcmyk")

T = 65
N_s = 22
A = 300
m = 7.64
R = 8.314
n_e = 2
F = 96485

r_1 = 5.95e-1
r_2 = -3.41e-3
r_3 = -1.07
r_4 = 27.08

s_1_a = 2.52e-2
s_2_a = -2.34e-4
s_3_a = 3.18e-6
t_1_a = 5.46e-2
t_2_a = -2.46e-3
t_3_a = 5.21e-5

s_1_c = 1.10e-1
s_2_c = -1.65e-3
s_3_c = 2.28e-5
t_1_c = 45.7
t_2_c = 7.78e-1
t_3_c = -1.06e-2

def celsius_to_kelvin(T_in_celsius):
    return T_in_celsius + 273.15

def calculate_vapour_pressure_of_water(T):
    T_kelvin = celsius_to_kelvin(T)
    return np.exp(81.6179 - (7699.68/T_kelvin) - 10.9*np.log(T_kelvin) + 9.5891*pow(10, -3)*T_kelvin)

def calculate_parameters_based_on_molal_concentration(m):
    a = -0.0151*m -(1.6788*pow(10, -3)*pow(m, 2)) + (2.2588*pow(10, -5)*pow(m, 3))
    b = 1 - 1.2062*pow(10, -3)*m + 5.6024*pow(10, -4)*pow(m, 2) - 7.8228*pow(10, -6)*pow(m, 3)
    return a, b

def calculate_saturated_vapour_pressure_of_water(T):
    balej = 35.4462 - 3343.93/T - 10.9*np.log10(T) + 4.1645e-3*T
    return pow(10, balej)

def calculate_vapour_pressure_electrolyte(T, m):
    a, b = calculate_parameters_based_on_molal_concentration(m)
    return (pow(10, a))*(pow(calculate_saturated_vapour_pressure_of_water(T), b))

def calculate_standard_rev_voltage(T):
    return 1.5184 - (1.5421*(10.0 ** (-3))*T) + (9.526*(10.0 ** (-5))*T*np.log(T)) + (9.84*(10.0 ** (-8))*T)

def calculate_rev_voltage(T, p, N_s, m, R, n_e, F):
    T_kelvin = celsius_to_kelvin(T)
    p_sv = calculate_saturated_vapour_pressure_of_water(T_kelvin)
    p_sv_el = calculate_vapour_pressure_electrolyte(T_kelvin, m)
    u_rev_0 = calculate_standard_rev_voltage(T_kelvin)
    f = R/(n_e*F)
    return N_s*(u_rev_0 + f*T_kelvin*np.log(((p - p_sv_el)**(3/2)*p_sv)/(p_sv_el)))

def calculate_ohmic_resistance(T, r_1, r_2, r_3, r_4, N_s, A):
    return N_s*((r_1 + r_2*T + (r_3/T) + (r_4/(T**2)))/A)

def calculate_parameter_s(T, s_1, s_2, s_3):
    return s_1 + s_2*T + s_3*pow(T, 2)

def calculate_parameter_t(T, t_1, t_2, t_3):
     return t_1 + t_2*T + t_3*pow(T, 2)

def calculate_activation_overvoltage(T, s_1, s_2, s_3, t_1, t_2, t_3, i_E, N_s):
    return N_s * calculate_parameter_s(T, s_1, s_2, s_3) * np.log((i_E/calculate_parameter_t(T, t_1, t_2, t_3)) + 1)

#temperature_list = [15, 25, 35, 45, 55, 65]
pressure_list = [15, 35, 65]
u_E = []
i_E = np.linspace(0, 120, 500)

for pressure in pressure_list:
    u_rev = calculate_rev_voltage(T, pressure, N_s, m, R, n_e, F)
    u_ohm = i_E * calculate_ohmic_resistance(T, r_1, r_2, r_3, r_4, N_s, A)
    u_act_c = calculate_activation_overvoltage(T, s_1_a, s_2_a, s_3_a, t_1_a, t_2_a, t_3_a, i_E, N_s)
    u_act_a = calculate_activation_overvoltage(T, s_1_c, s_2_c, s_3_c, t_1_c, t_2_c, t_3_c, i_E, N_s)
    u_E.append(u_rev + u_ohm + u_act_c +  u_act_a)


for index, t in enumerate(pressure_list):
    plt.plot(i_E, u_E[index], label=f'${t}$ bar', lw=2.5)

plt.xlabel('Current (A)')
plt.ylabel('Voltage (V)')
plt.xlim([0, 120])
plt.ylim([26, 42])
plt.xticks(np.arange(0, 121, 20))
plt.yticks(np.arange(26, 43, 2))
plt.title('Operating temperature: $65^\circ$ C')
plt.legend(loc='lower right')
plt.grid(True)

plt.show()