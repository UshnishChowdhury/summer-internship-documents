import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

T = 45
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


def calculate_parameter_s(T, s_1, s_2, s_3):
    return s_1 + s_2*T + s_3*pow(T, 2)

def calculate_parameter_t(T, t_1, t_2, t_3):
     return t_1 + t_2*T + t_3*pow(T, 2)

def calculate_activation_overvoltage(T, s_1, s_2, s_3, t_1, t_2, t_3, i_E, N_s):
    return N_s * calculate_parameter_s(T, s_1, s_2, s_3) * np.log((i_E/calculate_parameter_t(T, t_1, t_2, t_3)) + 1)

def calculate_activation_resistance(T, s_1, s_2, s_3, t_1, t_2, t_3, i_E, N_s):
    u_act = calculate_activation_overvoltage(T, s_1, s_2, s_3, t_1, t_2, t_3, i_E, N_s)
    return u_act/i_E

def calculate_ohmic_resistance(T, r_1, r_2, r_3, r_4, N_s, A):
    return N_s*((r_1 + r_2*T + (r_3/T) + (r_4/(T**2)))/A)

i_E_list = [10, 40, 80]

R_ohm_E = calculate_ohmic_resistance(T, r_1, r_2, r_3, r_4, N_s, A)
C_dl_a_E = 0.6375
C_dl_c_E = 0.0307

frequencies = np.logspace(-1, 3, num=100)
omega = 2 * np.pi * frequencies

Z_real = []
Z_imag = []

for i_E in i_E_list:
    R_act_a_E = calculate_activation_resistance(T, s_1_a, s_2_a, s_3_a, t_1_a, t_2_a, t_3_a, i_E, N_s)
    R_act_c_E = calculate_activation_resistance(T, s_1_c, s_2_c, s_3_c, t_1_c, t_2_c, t_3_c, i_E, N_s)
    print(R_act_a_E)
    Z_real.append(R_ohm_E + (R_act_a_E/(1 + np.pow(omega*R_act_a_E*C_dl_a_E, 2))) + (R_act_c_E/(1 + np.pow(omega*R_act_c_E*C_dl_c_E, 2))))
    Z_imag.append(-((omega*C_dl_a_E*np.pow(R_act_a_E, 2))/(1 + np.pow(omega*R_act_a_E*C_dl_a_E, 2)) 
                   + (omega*C_dl_c_E*np.pow(R_act_c_E, 2))/(1 + np.pow(omega*R_act_c_E*C_dl_c_E, 2))))

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# print(i_E_list)
# for index, i_E in enumerate(i_E_list):
#     ax.plot(i_E, Z_real[index], Z_imag[index])

# ax.set_xlabel('Current (A)')
# ax.set_ylabel('Real axis ($\Omega$)')
# ax.set_zlabel('Imaginary axis ($\Omega$)')

# plt.xticks(np.arange(20, 121, 40))


for index, i_E in enumerate(i_E_list):
    plt.plot(Z_real[index], Z_imag[index], label=f'{i_E} A', lw=2)

plt.xlabel('$Re\{ Z(j\omega) \} (\Omega)$')
plt.ylabel('$Im\{ Z(j\omega) \} (\Omega)$')
plt.xlim([0.03, 0.14])
plt.xticks(np.arange(0.03, 0.14, 0.01))
plt.ylim([-0.040, 0.000])
plt.legend(loc='lower right')
plt.grid(True)
plt.show()