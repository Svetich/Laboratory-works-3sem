from math import (sqrt, pi, log)

import numpy as np
import matplotlib.pyplot as plt
import os

import help_laba

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

B_I_file = open('Калибровка.txt')
B_I_str = B_I_file.read().splitlines()[1:]

B = np.zeros(len(B_I_str))
I = np.zeros(len(B_I_str))

for i in range(len(B_I_str)):
    data_I, data_B = B_I_str[i].split('\t')
    B[i] = float(data_B) * 10 ** (-3)
    I[i] = float(data_I)

b_BI, a_BI, sigma_BIb, epsilon_BIb = help_laba.LSM(I, B, len(I))

x_BI = np.linspace(0.04, 0.6, 100)
y_BI = b_BI * x_BI + a_BI

plt.plot(I, B, '.', markerfacecolor='crimson', markeredgecolor='crimson',
         markersize=10, markeredgewidth=1)
plt.plot(x_BI, y_BI, '-', color='dodgerblue', linewidth=1)
plt.xlabel(r'$I_{m}$, А', size=16)
plt.ylabel(r'$B$, Тл', size=16)
plt.show()
plt.close()

print('Коэф. пересчёта = ' + str(b_BI) + '+-' + str(sigma_BIb))

disk1_file = open('Диск1.txt')
disk1_str = disk1_file.read().splitlines()[1:]

Id = np.zeros(len(disk1_str))
Ud1 = np.zeros(len(disk1_str))

for i in range(len(disk1_str)):
    data_Id1, data_Ud1 = disk1_str[i].split('\t')
    Id[i] = float(data_Id1)
    Ud1[i] = float(data_Ud1)

disk2_file = open('Диск2.txt')
disk2_str = disk2_file.read().splitlines()[1:]

Ud2 = np.zeros(len(disk2_str))

for i in range(len(disk2_str)):
    data_Id2, data_Ud2 = disk2_str[i].split('\t')
    Ud2[i] = float(data_Ud2)

Ud = (Ud1 + Ud2) / 2

Id0 = 0.022
sigma_Ido = 0.001
epsilon_Ido = sigma_Ido / Id0

U0d = 0.29
sigma_Ud0 = 0.001
epsilon_Ud0 = sigma_Ud0 / U0d

Bd = b_BI * Id
Uyd = (Ud - U0d) / U0d

Bd_line = np.array([Bd[0], Bd[1], Bd[2], Bd[3]])
Uyd_line = np.array([Uyd[0], Uyd[1], Uyd[2], Uyd[3]])

b_d, a_d, sigmab_d, epsilonb_d = help_laba.LSM(Bd ** 2, Uyd, len(Uyd))

sigmaa_d = sigmab_d * sqrt(np.mean(Bd_line ** 4))
epsilona_d = sigmaa_d / a_d

x_d = np.linspace(0, 0.05, 100)
y_d = x_d * b_d + a_d

print('Коэф. наклона прямолинейного уч-ка гр-ка для диска Корбино = ' + str(b_d) + '+-0' + str(sigmab_d))

plate1_file = open('Пластина1.txt')
plate1_str = plate1_file.read().splitlines()[1:]

Ip1 = np.zeros(len(plate1_str))
Up1 = np.zeros(len(plate1_str))

for i in range(len(plate1_str)):
    data_Ip1, data_Up1 = plate1_str[i].split('\t')
    Ip1[i] = float(data_Ip1)
    Up1[i] = float(data_Up1)

Ip0 = 0.010
sigma_Ipo = 0.001
epsilon_Ipo = sigma_Ipo / Ip0

Up0 = 0.675
sigma_Up0 = 0.001
epsilon_Up0 = sigma_Up0 / Up0

Bp1 = b_BI * Ip1
Uyp1 = (Up1 - Up0) / Up0

plate2_file = open('Пластина2.txt')
plate2_str = plate2_file.read().splitlines()[1:]

Ip2 = np.zeros(len(plate2_str))
Up2 = np.zeros(len(plate2_str))

for i in range(len(plate2_str)):
    data_Ip2, data_Up2 = plate2_str[i].split('\t')
    Ip2[i] = float(data_Ip2)
    Up2[i] = float(data_Up2)

Bp2 = b_BI * Ip2
Uyp2 = (Up2 - Up0) / Up0


plt.plot(Bd ** 2, Uyd, '.', label='Диск Корбино', markerfacecolor='blue', markeredgecolor='blue',
         markersize=10, markeredgewidth=1)
plt.plot(Bp1 ** 2, Uyp1, '.', label='Пластина, средняя сторона вдоль поля', markerfacecolor='red', markeredgecolor='red',
         markersize=10, markeredgewidth=1)
plt.plot(Bp2 ** 2, Uyp2, '.', label='Пластина, средняя сторона поперёк поля', markerfacecolor='orange', markeredgecolor='orange',
         markersize=10, markeredgewidth=1)
plt.plot(x_d, y_d, '-', color='dodgerblue', linewidth=1)
plt.xlabel(r'$B^2$, мТл$^2$', size=16)
plt.ylabel(r'$\frac{U - U_{0}}{U_{0}}$', size=16)
plt.legend()
# plt.show()
plt.close()


mu = sqrt(b_d)
sigma_mu = mu * sqrt((1 / 2 * epsilonb_d) ** 2)
epsilon_mu = sigma_mu / mu

mu_th = 7.7
eps_mu = (mu_th - mu) / mu_th

print('Подвижность посителей = ' + str(mu) + '+-' + str(sigma_mu) + ' ' + str(eps_mu))

R0 = U0d / Id0
sigma_R0 = R0 * sqrt(epsilon_Ido ** 2 + epsilon_Ud0 ** 2)
epsilon_R0 = sigma_R0 / R0

print('Сопротивление диска в отсутствие магнитного поля = ' + str(R0) + '+-' + str(sigma_R0))

D = 18
sigmaD = 0.1
epsilonD = sigmaD / D

d = 3
sigmad = 0.1
epsilond = sigmad / d

h = 1.8 * 10 ** (-3)
sigmah = 0.1 * 10 ** (-3)
epsilonh = sigmah / h

rho0 = R0 * 2 * pi * h / log(D / d)
sigma_rho0 = rho0 * sqrt(epsilon_R0 ** 2 + epsilonh ** 2 + epsilonD ** 2 + epsilond ** 2)
epsilon_rho0 = sigma_rho0 / rho0

print('Удельное сопротивление = ' + str(rho0) + '+-' + str(sigma_rho0))

sigma0 = 1 / rho0
sigma_sigma0 = sigma0 * sqrt(epsilon_rho0 ** 2)
epsilon_sigma0 = sigma_sigma0 / sigma0

print('Удельная проводимость = ' + str(sigma0) + '+-' + str(sigma_sigma0))

n = sigma0 / (1.6 * 10 ** (-19) * mu)
sigma_n = n * sqrt(epsilon_sigma0 ** 2 + epsilon_mu ** 2)
epsilon_n = sigma_n / n

print('Концентрация носителей тока = ' + str(n) + '+-' + str(sigma_n))




