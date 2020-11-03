from math import sqrt
from math import pi

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


def mean_square_error(data):
    data_mean = np.mean(data)
    delta = data - data_mean
    sigma = sqrt(len(data) / (len(data) - 1) * sum(delta ** 2))
    return data_mean, sigma


p8 = open('p.8', 'r')
p8_str = p8.read().splitlines()[1:]

C = np.zeros(len(p8_str))
f = np.zeros(len(p8_str))
Uc = np.zeros(len(p8_str))
E = np.zeros(len(p8_str))
omega0 = np.zeros(len(p8_str))
L = np.zeros(len(p8_str))
Q = np.zeros(len(p8_str))
ro = np.zeros(len(p8_str))
Rsum = np.zeros(len(p8_str))
Rsmax = np.zeros(len(p8_str))
Rl = np.zeros(len(p8_str))
I = np.zeros(len(p8_str))

for i in range(len(p8_str)):
    data_C, data_f, data_Uc, data_E = p8_str[i].split('\t')
    C[i] = 10 ** (-9) * float(data_C.replace(',', '.'))
    f[i] = 10 ** 3 * float(data_f.replace(',', '.'))
    Uc[i] = float(data_Uc.replace(',', '.'))
    E[i] = float(data_E.replace(',', '.'))
    omega0[i] = 2 * pi * f[i]
    L[i] = 1 / (omega0[i] ** 2 * C[i])
    Q[i] = Uc[i] / E[i]
    ro[i] = sqrt(L[i] / C[i])
    Rsum[i] = ro[i] / Q[i]
    Rsmax[i] = 10 ** (-3) * ro[i]
    Rl[i] = Rsum[i] - Rsmax[i] - 3.45
    I[i] = E[i] / Rsum[i]


print(Rl)

L_mean, sigmaL = mean_square_error(L)
Rl_mean, sigmaRl = mean_square_error(Rl)


print('Среднее значение индуктивности = ' + str(L_mean) + '+-' + str(sigmaL))
print('Среднее значение активного сопротивления = ' + str(Rl_mean) + '+-' + str(sigmaRl))

p9 = open('АЧХ', 'r')
p9_str = p9.read().splitlines()[2:]

f1 = np.zeros(len(p9_str))  # С = 33,2нФ
Uc1 = np.zeros(len(p9_str))
f2 = np.zeros(len(p9_str))  # С = 57,2нФ
Uc2 = np.zeros(len(p9_str))

f1_f0 = np.zeros(len(p9_str))
u1_u0 = np.zeros(len(p9_str))
f2_f0 = np.zeros(len(p9_str))
u2_u0 = np.zeros(len(p9_str))

for i in range(len(p9_str)):
    data_f1, data_u1, data_f2, data_u2 = p9_str[i].split('\t')
    f1[i] = float(data_f1.replace(',', '.'))
    Uc1[i] = float(data_u1.replace(',', '.'))
    f2[i] = float(data_f2.replace(',', '.'))
    Uc2[i] = float(data_u2.replace(',', '.'))
    f1_f0[i] = f1[i] / 28.5
    f2_f0[i] = f2[i] / 21.5
    u1_u0[i] = Uc1[i] / 5.07
    u2_u0[i] = Uc2[i] / 5.7


plt.plot(f1, Uc1, '.', color='orange')
plt.plot(f2, Uc2, '.', color='blue')
plt.show()
plt.close()

plt.plot(f1_f0, u1_u0, '.', color='orange')
plt.plot(f2_f0, u2_u0, '.', color='blue')
plt.show()
plt.close()









