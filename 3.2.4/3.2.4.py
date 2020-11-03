from math import (pi, sqrt, log)

import numpy as np
import matplotlib.pyplot as plt

import help_laba

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


def error_ln(x, y, sigma_x, sigma_y):
    sigma_ln = sqrt((sigma_x / x) ** 2 + (sigma_y / y) ** 2)
    return sigma_ln


T0 = 0.01
x0 = 10
L = 200

Cnx_file = open('2. Период свободных колебаний')
Cnx_str = Cnx_file.read().splitlines()[1:]

C = np.zeros(len(Cnx_str))
n = np.zeros(len(Cnx_str))
x = np.zeros(len(Cnx_str))
Tth = np.zeros(len(Cnx_str))

for i in range(len(Cnx_str)):
    dataC, dataT0, datax0, datan, datax = Cnx_str[i].split('\t')
    C[i] = float(dataC)
    n[i] = float(datan)
    x[i] = float(datax)
    Tth[i] = 2 * pi * sqrt(L * 10 ** (-9) * C[i])

Texp = T0 * x / (n * x0)

b_tt, a_tt, sigma_tt, epsilon_tt = help_laba.LSM(Tth, Texp, len(Tth))

print('Коэф. наклона прямой Texp(Tth) = ' + str(b_tt) + '+-' + str(sigma_tt))

x_tt = np.linspace(0.0003, 0.003, 100)
y_tt = b_tt * x_tt + a_tt

plt.plot(x_tt, y_tt, '-', label='1', color='royalblue', linewidth=1)
plt.plot(Tth, Texp, '.', label='2', markerfacecolor='crimson', markeredgecolor='crimson',
         markersize=10, markeredgewidth=1)
plt.legend()
plt.xlabel(r'$T_{th}$, c', size=18)
plt.ylabel(r'$T_{exp}$, c', size=18)
# plt.show()
plt.close()


Rl = 38.957 * 10 ** (-3)

RnU_file = open('3. Декремент затухания')
RnU_str = RnU_file.read().splitlines()[1:]

R = np.zeros(len(RnU_str))
n3 = np.zeros(len(RnU_str))
Uk = np.zeros(len(RnU_str))
Ukn = np.zeros(len(RnU_str))
tetta1 = np.zeros(len(RnU_str))

for i in range(len(RnU_str)):
    dataR, datan3, dataUk, dataUkn = RnU_str[i].split('\t')
    R[i] = float(dataR)
    n3[i] = float(datan3)
    Uk[i] = float(dataUk)
    Ukn[i] = float(dataUkn)
    tetta1[i] = 1 / n3[i] * log(Uk[i] / Ukn[i])

Rsum = R + Rl

xexp_cr = 1 / Rsum ** 2
yexp_cr = 1 / tetta1 ** 2

Rcrexp, aRcr, sigma_Rcrexp, epsilon_Rcrexp = help_laba.LSM(xexp_cr, yexp_cr, len(xexp_cr))

print('Экспериментальное Rcr = ' + str(sqrt(Rcrexp) * 2 * pi) + '+-' + str(Rcrexp * (1/2 * epsilon_Rcrexp)))

Q_min = 2 * pi / max(tetta1)
Q_max = 2 * pi / min(tetta1)

sigma_Q_min = Q_min * sqrt((1 / n[5]) ** 2 + (error_ln(Uk[5], Ukn[5], 0.1, 0.1)/log(Uk[5] / Ukn[5])) ** 2)
sigma_Q_max = Q_max * sqrt((1 / n[0]) ** 2 + (error_ln(Uk[0], Ukn[0], 0.1, 0.1)/log(Uk[0] / Ukn[0])) ** 2)

epsilon_Q_min = sigma_Q_min / Q_min
epsilon_Q_max = sigma_Q_max / Q_max

print('Минимальная добротность по графику = ' + str(Q_min) + '+-' + str(sigma_Q_min))
print('Максимальная добротность по графику = ' + str(Q_max) + '+-' + str(sigma_Q_max))

xexp_cr_line = np.linspace(0.7, 0.05, 100)
yexp_cr_line = xexp_cr_line * Rcrexp + aRcr

plt.plot(xexp_cr_line, yexp_cr_line, '-', label='1', color='royalblue', linewidth=1)
plt.plot(xexp_cr, yexp_cr, '.', label='2', markerfacecolor='crimson', markeredgecolor='crimson',
         markersize=10, markeredgewidth=1)
plt.legend()
plt.xlabel(r'$1/R_{sum}^2$, 1/Ом$^2$', size=18)
plt.ylabel(r'$1/\Theta^2$', size=18)
plt.show()
plt.close()

Spiral_file = open('4. Фазовая плоскость')
Spiral_str = Spiral_file.read().splitlines()[1:]

R_sp = np.zeros(len(Spiral_str))
n_sp = np.zeros(len(Spiral_str))
xk = np.zeros(len(Spiral_str))
xkn = np.zeros(len(Spiral_str))
tetta2 = np.zeros(len(Spiral_str))

for i in range(len(Spiral_str)):
    dataR_sp, datan_sp, dataxk, dataxkn = Spiral_str[i].split('\t')
    R_sp[i] = float(dataR_sp)
    n_sp[i] = float(datan_sp)
    xk[i] = float(dataxk)
    xkn[i] = float(dataxkn)
    tetta2[i] = 1 / n_sp[i] * log(xkn[i] / xk[i])

Q_sp_min = pi / max(tetta2)
Q_sp_max = pi / min(tetta2)

sigma_Q_min_sp = Q_sp_min * sqrt((1 / n_sp[5]) ** 2 + (error_ln(xk[5], xkn[5], 0.1, 0.1)/log(xk[5] / xkn[5])) ** 2)
sigma_Q_max_sp = Q_sp_max * sqrt((1 / n_sp[0]) ** 2 + (error_ln(xk[0], xkn[0], 0.1, 0.1)/log(xk[0] / xkn[0])) ** 2)

epsilon_Q_min_sp = sigma_Q_min_sp / Q_sp_min
epsilon_Q_max_sp = sigma_Q_max_sp / Q_sp_max

print('Минимальная добротность по спирали = ' + str(Q_sp_min) + '+-' + str(sigma_Q_min_sp))
print('Максимальная добротность по спирали = ' + str(Q_sp_max) + '+-' + str(sigma_Q_max_sp))

R_cr_th = 2 * sqrt(200.6 * 10 ** (-3) / (0.005 * 10 ** (-6)))
sigma_Rcrth = R_cr_th * sqrt((1/2 * (0.01 / 200.6)) ** 2 + (1/2 * 0.001 / 0.005) ** 2)

print('Теоретическое значение Rcrth = ' + str(R_cr_th / 1000) + '+-' + str(sigma_Rcrth / 1000))

Q_max_th = 1 / 2 * sqrt((R_cr_th / (1000 * R[0])) ** 2 - 1)
Q_min_th = 1 / 2 * sqrt((R_cr_th / (1000 * R[5])) ** 2 - 1)

sigma_Q_max_th = Q_max_th * sqrt((sigma_Rcrth / R_cr_th) ** 2 + (0.1 / R[0]) ** 2)
sigma_Q_min_th = Q_min_th * sqrt((sigma_Rcrth / R_cr_th) ** 2 + (0.1 / R[5]) ** 2)

print('Максимальная теоретическая добротность = ' + str(Q_max_th) + '+-' + str(sigma_Q_max_th))
print('Минимальная теоретическая добротность = ' + str(Q_min_th) + '+-' + str(sigma_Q_min_th))