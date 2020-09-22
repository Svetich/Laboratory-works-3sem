from math import sqrt
from math import pi
from math import atan

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


B_p_mes = 114 * 10  # мТл
sigma_B_p_mes = 0.001  # мТл

m = np.array([0.495, 0.851])  # массы каждого шарика, г
m_mean = np.mean(m)  # средняя масса шариков
sigma_m = 10 ** (-3)
epsil_m = sigma_m / m_mean

print('Средняя масса шариков = ' + str(m_mean) + ' г')

d = np.array([0.5, 0.568])  # диаметры шариков, см
d_mean = np.mean(d)
sigma_d_1 = 10 ** (-3)
sigma_d = 2 * sigma_d_1
epsil_d = sigma_d / d_mean

print('Средний диаметр шариков = ' + str(d_mean) + ' см')

r_max_mes = 1.557  # см
r_max = r_max_mes + d_mean
sigma_r_max = 10 ** (-2)

print('r_max = ' + str(r_max) + ' см')

g = 9.8 * 100  # см/с^2 (Гал)

magnetic_moment = sqrt(m_mean * g * (r_max ** 4) / 6)  # эрг/Гс
sigma_magnetic_moment = magnetic_moment * sqrt(
    ((0.5 * epsil_m) ** 2) + ((2 * (sigma_r_max / r_max)) ** 2))  # магнитный момент
epsil_magnetic_moment = sigma_magnetic_moment/magnetic_moment

print('Магнитный момент = ' + str(magnetic_moment) + '+-' +
      str(sigma_magnetic_moment) + ' эрг/Гс, ' +
      str(epsil_magnetic_moment))


V1 = 4/3 * pi * ((d[0] / 2) ** 3)  # объём 1го шарика, см^3
V2 = 4/3 * pi * ((d[1] / 2) ** 3)  # объём 2го шарика, см^3
V_all = [V1, V2]
V_mean = np.mean(V_all)  # средний объём
sigma_V1 = V1 * sqrt((3 * sigma_d/d[0]) ** 2)
sigma_V2 = V2 * sqrt((3 * sigma_d/d[1]) ** 2)
sigma_V = sigma_V1 + sigma_V2
epsil_V = sigma_V/V_mean

M = magnetic_moment / V_mean  # намагниченность, эрг/(Гс * см^3)
sigma_M = M * sqrt(epsil_magnetic_moment ** 2 + epsil_V ** 2)
epsil_M = sigma_M / M

print('Намагниченность = ' + str(M) + '+-' + str(sigma_M) + ' эрг/(Гс * см^3), ' + str(epsil_M))

B_r = 4 * pi * M
sigma_B_r = B_r * epsil_M
epsil_B_r = sigma_B_r/B_r

print('Остаточная индукция = ' + str(B_r) + '+-' + str(sigma_B_r) + 'Гс, ' + str(epsil_B_r))

teor_B_r = 1.22 * 10 ** 4
error_B_r = (teor_B_r - B_r)/teor_B_r

print('Расхождение Br = ' + str(error_B_r * 100))

B_p = 2/3 * B_r
sigma_B_p = B_p * epsil_B_r
epsil_B_p = sigma_B_p / B_p

print('Индукция у полюсов = ' + str(B_p) + '+-' + str(sigma_B_p) + 'Гс, ' + str(epsil_B_p))

error_B_p = (B_p_mes - B_p)/B_p_mes

print('Расхождения Bp = ' + str(error_B_p * 100))

T_line = 30.20 / 10
sigma_T_line = T_line * 0.01 / 30.20
T_ring = 55.31 / 10
sigma_T_ring = T_ring * 0.01 / 55.31

print('Период колебаний линии = ' + str(T_line) + '+-' + str(sigma_T_line))
print('Период колебаний кольца = ' + str(T_ring) + '+-' + str(sigma_T_ring))

t = np.array([33.28, 29.25, 28.23, 26.21, 22.18, 19.15, 16.13, 14.10, 12.08, 10.06])
N = np.linspace(10, 10, 10)
n = np.array([12, 11, 10, 9, 8, 7, 6, 5, 4, 3])
T = t / N
b_Tn, sigma_b_Tn, epsilon_b_Tn = LSM(n, T, 10)

print('Коэф. наклона T(n) = ' + str(b_Tn) + '+-' + str(sigma_b_Tn) +
      ', ' + str(epsilon_b_Tn))

B_hor = (m_mean * (d_mean * pi) ** 2)/\
        (3 * magnetic_moment * (b_Tn ** 2))
sigma_B_hor = B_hor * sqrt(epsil_m ** 2 + (2 * epsil_d) ** 2
                           + (2 * epsilon_b_Tn) ** 2)
epsil_B_hor = sigma_B_hor/B_hor

print('Горизонтальная сост-я = ' + str(B_hor) + '+-' + str(sigma_B_hor) +
      ' Гс, ' + str(epsil_B_hor))

x_Tn = np.linspace(0, 15, 15)
y_Tn = b_Tn * x_Tn

plt.plot(n, T, '.', label='экспериментальные точки', markerfacecolor='red', markeredgecolor='red',
         markersize=10, markeredgewidth=1)
plt.plot(x_Tn, y_Tn, label='аппрокс. точки', color='black', linewidth=1)
plt.xlabel(r'$n$', size=18)
plt.ylabel(r'$T$, c', size=18)
# plt.legend()
# plt.show()
# plt.savefig('T(n).png')
plt.close()

n_even = np.array([8, 6, 4])
m_load = np.array([0.314, 0.221, 0.172])
g_array = np.linspace(9.8, 9.8, 3)
power = m_load * g_array * 100
shoulders = np.array([d_mean * 3, d_mean * 2, d_mean])
moment_power = shoulders * power

n_bad = np.array([12, 10])  # плохие точки
m_bad = np.array([0.169, 0.172])
g_bad = np.linspace(9.8, 9.8, 2)
power_bad = m_bad * g_bad
shoulders_bad = np.array([d_mean * 5, d_mean * 4])
moment_power_bad = shoulders_bad * power_bad

b_Mn, sigma_b_Mn, epsilon_b_Mn = LSM(n_even, moment_power, 3)

print('Коэф. наклона M(n) = ' + str(b_Mn) + '+-' + str(sigma_b_Mn) +
      ', ' + str(epsilon_b_Mn))

B_ver = b_Mn / magnetic_moment
sigma_B_ver = B_ver * sqrt(epsilon_b_Mn ** 2 + epsil_magnetic_moment ** 2)
epsil_B_ver = sigma_B_ver / B_ver

print('Вертикальная сост-я = ' + str(B_ver) + '+-' + str(sigma_B_ver) +
      ' Гс, ' + str(epsil_B_ver))

tg_beta = B_ver/B_hor
sigma_tg_beta = tg_beta * sqrt(epsil_B_ver**2 + epsil_B_hor**2)
epsil_tg_beta = sigma_tg_beta / tg_beta

print('Магнитное наклонение = ' + str(atan(tg_beta)) + ' (' + str(tg_beta) + '+-' + str(sigma_tg_beta) + str(epsil_tg_beta) + ')')

B = sqrt(B_ver ** 2 + B_hor ** 2)
sigma_B = B * sqrt((0.5 * epsil_B_hor) ** 2 + (0.5 * epsil_B_ver) ** 2)
epsil_B = sigma_B/B

print('Индукция МП Земли = ' + str(B) + '+-' + str(sigma_B) + 'Гс, ' +str(epsil_B))

x_Mn = np.linspace(0, 10, 100)
y_Mn = b_Mn * x_Mn

plt.plot(n_even, moment_power, '.', label='экспериментальные точки', markerfacecolor='red', markeredgecolor='red',
         markersize=10, markeredgewidth=1)
plt.plot(n_bad, moment_power_bad, '.', label='выпавшие точки', markerfacecolor='blue', markeredgecolor='blue',
         markersize=10, markeredgewidth=1)
plt.plot(x_Mn, y_Mn, label='аппрокс. точки', color='black', linewidth=1)
plt.xlabel(r'$n$', size=18)
plt.ylabel(r'$M$, $дин \cdot см$', size=18)
# plt.savefig('M(n).png')
# plt.show()
plt.close()


