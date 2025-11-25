import matplotlib.pyplot as plt
import numpy as np

N_s = 22
A = 300

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


def calculate_ohmic_resistance(T, r_1, r_2, r_3, r_4, N_s, A):
    return N_s*((r_1 + r_2*T + (r_3/T) + (r_4/(T**2)))/A)

def celsius_to_kelvin(T_in_celsius):
    return T_in_celsius + 273.15

def calculate_vapour_pressure_of_water(T):
    T_kelvin = celsius_to_kelvin(T)
    return np.exp(81.6179 - (7699.68/T_kelvin) - 10.9*np.log(T_kelvin) + 9.5891*pow(10, -3)*T_kelvin)

def calculate_parameters_based_on_molal_concentration(m):
    a = -0.0151*m -(1.6788*pow(10, -3)*pow(m, 2)) + (2.2588*pow(10, -5)*pow(m, 3))
    b = 1 - 1.2062*pow(10, -3)*m + 5.6024*pow(10, -4)*pow(m, 2) - 7.8228*pow(10, -6)*pow(m, 3)
    return a, b

def calculate_vapour_pressure_KOH_solution(a, b, p_v_h2o):
    return np.exp(2.302*a + b*np.log(p_v_h2o))

def calculate_water_activity_in_KOH_solution(m, T):
    T_kelvin = celsius_to_kelvin(T)
    return np.exp(-0.05192*m + 0.003302*pow(m, 2) + (3.177*m - 2.131*pow(m, 2))/T_kelvin)

def calculate_standard_rev_voltage(T):
    T_kelvin = celsius_to_kelvin(T)
    return 1.5184 - (1.5421*(10.0 ** (-3))*T_kelvin) + (9.526*(10.0 ** (-5))*T_kelvin*np.log(T_kelvin)) + (9.84*(10.0 ** (-8))*T_kelvin)

def calculate_rev_voltage(T, p, N_s, m, R, n_e, F):
    T_kelvin = celsius_to_kelvin(T)
    a, b = calculate_parameters_based_on_molal_concentration(m)
    p_v_h2o = calculate_vapour_pressure_of_water(T)
    pressure_error = p - calculate_vapour_pressure_KOH_solution(a, b, p_v_h2o)
    return N_s*(calculate_standard_rev_voltage(T) +
           R*T_kelvin*np.log(pressure_error*np.sqrt(pressure_error)/calculate_water_activity_in_KOH_solution(m, T))/(n_e*F))


def calculate_parameter_s(T, s_1, s_2, s_3):
    return s_1 + s_2*T + s_3*pow(T, 2)

def calculate_parameter_t(T, t_1, t_2, t_3):
     return t_1 + t_2*T + t_3*pow(T, 2)

def calculate_activation_overvoltage(T, s_1, s_2, s_3, t_1, t_2, t_3, i_E, N_s):
    return N_s * calculate_parameter_s(T, s_1, s_2, s_3) * np.log((i_E/calculate_parameter_t(T, t_1, t_2, t_3)) + 1)

#temperature_list = [15, 25, 35, 45, 55, 65]
temperature_list = [15, 35, 65]
u = []
i_E = np.linspace(0, 120, 500)

for temperature in temperature_list:
    u_rev = calculate_rev_voltage(temperature, 25, 22, 7.64, 8.314, 2, 96485)
    u_ohm = i_E * calculate_ohmic_resistance(temperature, 5.95*pow(10,-1), -3.41*pow(10,-3), -1.07, 27.08, 22, 300)
    u_act_c = calculate_activation_overvoltage(temperature, 1.10e-1, -1.65e-3, 2.28e-5, 45.7, 7.78e-1, -1.06e-2, i_E, 22)
    u_act_a = calculate_activation_overvoltage(temperature, 2.52e-2, -2.34e-4, 3.18e-6, 5.46e-2, -2.46e-3, 5.21e-5, i_E, 22)
    u.append(u_rev + u_ohm + u_act_c +  u_act_a)

plt.plot(i_E, u[0], label="$15^\circ$ C", lw=2)

#plt.plot(i_E, u[1], label="$25^\circ$ C")

plt.plot(i_E, u[1], label="$35^\circ$ C", lw=2)

#plt.plot(i_E, u[3], label="$45^\circ$ C")

#plt.plot(i_E, u[4], label="$55^\circ$ C")

plt.plot(i_E, u[2], label="$65^\circ$ C", lw=2)

plt.xlabel('Current (A)')
plt.ylabel('Voltage (V)')

plt.xticks(np.arange(0, 121, 20))
plt.yticks(np.arange(26, 43, 2))

plt.legend(loc='lower right')
plt.grid(True)

plt.show()