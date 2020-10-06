from math import sqrt
from math import pi
from math import log

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


def LSM(x, y, n):
    x_y = np.mean(x * y)
    x_ = np.mean(x)
    y_ = np.mean(y)
    x_2 = np.mean(x ** 2)
    y_2 = np.mean(y ** 2)
    b = (x_y - (x_ * y_))/\
        (x_2 - (x_ ** 2))
    sigma_b = 1/sqrt(n) * sqrt((y_2 - (y_ ** 2))/
                               (x_2 - (x_ ** 2)) - b ** 2)
    epsilon = sigma_b / b
    return b, sigma_b, epsilon


VAHfileP = open(r'Характеристика стабилитрона (на повышение).txt')
VAHp = VAHfileP.read().splitlines()[1:]

Up = np.zeros(len(VAHp))
Ip = np.zeros(len(VAHp))

for i in range(len(VAHp)):
    data_Up1, data_Ip1 = VAHp[i].split('\t\t')
    data_Up2 = list(map(float, data_Up1.split(', ')))
    data_Ip2 = list(map(float, data_Ip1.split(', ')))
    Up[i] = np.mean(data_Up2)
    Ip[i] = np.mean(data_Ip2)


VAHfileM = open(r'Характеристика стабилитрона (на понижение).txt')
VAHm = VAHfileM.read().splitlines()[1:]

Um = np.zeros(len(VAHm))
Im = np.zeros(len(VAHm))

for i in range(len(VAHm)):
    data_Um1, data_Im1 = VAHm[i].split('\t\t')
    data_Um2 = list(map(float, data_Um1.split(', ')))
    data_Im2 = list(map(float, data_Im1.split(', ')))
    Um[i] = np.mean(data_Um2)
    Im[i] = np.mean(data_Im2)

VAHfileall = open(r'Характеристика стабилитрона (вместе).txt')
VAHall = VAHfileall.read().splitlines()[1:]

Uall = np.zeros(len(VAHall))
Iall = np.zeros(len(VAHall))

for i in range(len(VAHall)):
    data_Uall1, data_Iall1 = VAHall[i].split('\t\t')
    data_Uall2 = list(map(float, data_Uall1.split(', ')))
    data_Iall2 = list(map(float, data_Iall1.split(', ')))
    Uall[i] = np.mean(data_Uall2)
    Iall[i] = np.mean(data_Iall2)

b_UI, sigma_bUI, epsilon_bUI = LSM(Uall, Iall, len(Iall))
a = (np.mean(Iall) - b_UI * np.mean(Uall))

print('Коэф. наклона прямой экспер. = ' + str(b_UI) + '+-' + str(sigma_bUI))

xU = np.linspace(79, 175, 150)
yU = xU * b_UI + a

plt.plot(xU, yU, '-', label='1', color='orange', linewidth=1)
plt.plot(Up, Ip, '.', label='2', markerfacecolor='red', markeredgecolor='red',
         markersize=10, markeredgewidth=1)
plt.plot(Um, Im, '.', label=3, markerfacecolor='blue', markeredgecolor='blue',
         markersize=10, markeredgewidth=1)
plt.legend()
plt.xlabel(r'$U$, B', size=18)
plt.ylabel(r'$I$, мА', size=18)
# plt.show()
# plt.savefig('UI.png')
plt.close()

r = 5.1
Up_r = Up - r * Ip
Um_r = Um - r * Im
Uall_r = Uall - r * Iall

plt.plot(Up_r, Ip, '.', label='1', markerfacecolor='red', markeredgecolor='red',
         markersize=10, markeredgewidth=1)
plt.plot(Um_r, Im, '.', label='2', markerfacecolor='blue', markeredgecolor='blue',
         markersize=10, markeredgewidth=1)
plt.legend()
plt.xlabel(r'$U$, B', size=18)
plt.ylabel(r'$I$, мА', size=18)
# plt.show()
# plt.savefig('UI.png')
plt.close()

U = 108.15
sigma_U = 1
epsilon_U = sigma_U / U

V1 = 83.51
sigma_system_V = 0.01
sigma_random_1 = sqrt(1/2 * ((83.43 - V1) ** 2 + (83.52 - V1) ** 2 + (83.38 - V1) ** 2))
sigma_V1 = sqrt(sigma_system_V ** 2 + sigma_random_1 ** 2)
epsilon_V1 = sigma_V1 / V1

print('Напряжение зажигания = ' + str(V1) + '+-' + str(sigma_V1) + ' (' + str(epsilon_V1) + ')')

V2 = 80.57
sigma_random_2 = sqrt(1/2 * ((80.44 - V2) ** 2 + (80.62 - V2) ** 2 + (80.65 - V2) ** 2))
sigma_V2 = sqrt(sigma_system_V ** 2 + sigma_random_2 ** 2)
epsilon_V2 = sigma_V2 / V2

print('Напряжение гашения = ' + str(V2) + '+-' + str(sigma_V2) + ' (' + str(epsilon_V2) + ')')

I1 = 2.196
sigma_system_I = 0.001
sigma_random_I1 = sqrt(1/2 * ((2.187 - I1) ** 2 + (2.2 - I1) ** 2 + (2.202 - I1) ** 2))
sigma_I1 = sqrt(sigma_system_I ** 2 + sigma_random_I1 ** 2)

print('Сила тока I1 = ' + str(I1) + '+-' + str(sigma_I1))

I2 = 1.684
sigma_random_I2 = sqrt(1/2 * ((1.687 - I2) ** 2 + (1.684 - I2) ** 2 + (1.682 - I2) ** 2))
sigma_I2 = sqrt(sigma_system_I ** 2 + sigma_random_I2 ** 2)

print('Сила тока I2 = ' + str(I2) + '+-' + str(sigma_I2))

const = log((108.15 - 80.57) / (100 - 83.51))
sigma_const = sqrt((sigma_V2 / (U - V2)) ** 2 + (sigma_V1 / (U - V1)) ** 2 +
                   (((V2 - V1) * sigma_U)/((U - V2) * (U - V1))) ** 2)
epsilon_const = sigma_const / const

R_one = 900
sigma_R_one = 1
epsilon_R_one = sigma_R_one / R_one

C_one = 0.05
sigma_C_one = 0.01
epsilon_C_one = sigma_C_one / C_one

nuC_file = open('Частота и ёмкость.txt')
nuC = nuC_file.read().splitlines()[1:]

nu_exp_C = np.zeros(len(nuC))
C = np.zeros(len(nuC))

for i in range(len(nuC)):
    data_C, data_nu = nuC[i].split('\t\t')
    nu_exp_C[i] = float(data_nu)
    C[i] = float(data_C)

T_exp_C = 1 / nu_exp_C

b_TexpC, sigma_TexpC, epsilon_TexpC = LSM(C, T_exp_C, len(C))
a_TexpC = (np.mean(T_exp_C) - b_TexpC * np.mean(C))

x_TexpC = np.linspace(2.5e-2, 5.5e-2, 100)
y_TexpC = x_TexpC * b_TexpC + a_TexpC

x_TthC = np.linspace(2.5e-2, 5.5e-2, 100)
y_TthC = 600 * (10 ** (-3)) * x_TthC * const

sigma_b_thC = 600 * (10 ** (-3)) * const * sqrt(epsilon_const ** 2 + epsilon_R_one ** 2)

plt.plot(x_TexpC, y_TexpC, '-', label='1', color='blue', linewidth=1)
plt.plot(5e-2, 1/28, '.', label='2', markerfacecolor='black', markeredgecolor='black',
         markersize=10, markeredgewidth=1)
plt.plot(C, T_exp_C, '.', label='3', markerfacecolor='red', markeredgecolor='red',
         markersize=10, markeredgewidth=1)
plt.plot(x_TthC, y_TthC, '-', label=4, color='orange', linewidth=1)
plt.legend()
plt.xlabel(r'$C$, мкФ', size=18)
plt.ylabel(r'$T_{exp}$, с', size=18)
# plt.show()
plt.close()

print('Коэф. наклона прямой Tэкс(C) = ' + str(b_TexpC) + '+-' + str(sigma_TexpC))
print('Коэф. наклона прямой Tтеор(С) = ' + str(600 * (10 ** (-3)) * const) + '+-' + str(sigma_b_thC))

nuR_file = open('Частота и сопротивление.txt')
nuR = nuR_file.read().splitlines()[1:]

R = np.zeros(len(nuR))
nu_exp_R = np.zeros(len(nuR))

for i in range(len(nuR)):
    data_R, data_nu1 = nuR[i].split('\t\t')
    R[i] = float(data_R)
    nu_exp_R[i] = float(data_nu1)

T_exp_R = 1 / nu_exp_R

b_TexpR, sigma_TexpR, epsilon_TexpR = LSM(R, T_exp_R, len(R))
a_TexpR = (np.mean(T_exp_R) - b_TexpR * np.mean(R))

x_TexpR = np.linspace(3.5e2, 9e2, 100)
y_TexpR = b_TexpR * x_TexpR + a_TexpR

x_TthR = np.linspace(3.5e2, 9e2, 100)
y_TthR = 0.05 * (10 ** (-3)) * x_TthR * const

sigma_b_thR = 0.05 * (10 ** (-3)) * const * sqrt(epsilon_const ** 2 + epsilon_C_one ** 2)

plt.plot(x_TexpR, y_TexpR, '-', label='1', color='blue', linewidth=1)
plt.plot(R, T_exp_R, '.', label='2', markerfacecolor='red', markeredgecolor='red',
         markersize=10, markeredgewidth=1)
plt.plot(x_TthR, y_TthR, '-', label='3', color='orange', linewidth=1)
plt.legend()
plt.xlabel(r'$R$, кОм', size=18)
plt.ylabel(r'$T_{exp}$, с', size=18)
# plt.show()
plt.close()

print('Коэф. наклона прямой Tэкс(R) = ' + str(b_TexpR) + '+-' + str(sigma_TexpR))
print('Коэф. наклона прямой Tтеор(R) = ' + str(0.05 * (10 ** (-3)) * const) + '+-' + str(sigma_b_thR))

Rcr_th = (U - V2) / I2
sigma_Rcr_th = Rcr_th * sqrt((sqrt(sigma_U ** 2 + sigma_V2 ** 2)/(U - V2)) ** 2 + (sigma_I2 / I2) ** 2)

print('Теоретическое критическое сопротивление = ' + str(Rcr_th) + '+-' + str(sigma_Rcr_th))


