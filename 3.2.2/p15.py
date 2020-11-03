import numpy as np
import matplotlib.pyplot as plt
import math

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


def approx(X, Y, d):
    m = len(X)
    xav = 0
    yav = 0
    x2av = 0
    y2av = 0
    xyav = 0
    # 1 approximation : y = kx + b
    if (d == 1):
        for n in range(len(X)):
            xyav = xyav + X[n] * Y[n] / m
            y2av = y2av + Y[n] * Y[n] / m
            x2av = x2av + X[n] * X[n] / m
            xav = xav + X[n] / m
            yav = yav + Y[n] / m
        k = ((xyav) - xav * yav) / (x2av - xav * xav)
        b = yav - k * xav
        d_k = 1 / math.sqrt(m) * math.sqrt(abs((y2av - yav * yav) / (x2av - xav * xav) - k ** 2))
        d_b = d_k * math.sqrt(x2av - xav * xav)
        final_data = [k, d_k, b, d_b]
    else:
        # 2 approximation : y = kx
        if (d == 2):
            for n in range(len(X)):
                xyav = xyav + X[n] * Y[n] / m
                y2av = y2av + Y[n] * Y[n] / m
                x2av = x2av + X[n] * X[n] / m
            k = (xyav) / (x2av)
            d_k = 1 / math.sqrt(m) * math.sqrt(abs((y2av) / (x2av) - k ** 2))
            final_data = [k, d_k]
    return final_data


f = [32.5, 28.5, 23.9, 21.5, 19.83, 17.96, 16.42]
err_f = [0 for i in range(len(f))]
for i in range(len(f)):
    err_f[i] = 0.03 * f[i]

RL = [6.51, 3.03, 2.89, 0.99, 1.01, 0.90, 1.11]
err_RL = [0 for i in range(len(RL))]
for i in range(len(RL)):
    err_RL[i] = 0.03 * RL[i]

plt.scatter(f, RL, s=5, color='navy')

final_data = approx(f, RL, 1)
k = final_data[0]
b = final_data[2]

t = np.arange(16, 36, 0.1)
Rl_mean = np.linspace(2.35, 2.35, 200)
s = [0 for i in range(len(t))]
for i in range(len(s)):
    s[i] = k * t[i] + b
plt.plot(t, Rl_mean, '--', color='black')
plt.plot(t, s, lw=1.0, color='#36B9DA')


plt.xlabel(r'$\nu_0$, кГц', size=15)
plt.ylabel(r'$R_{L}$, Ом', size=15)
kwargs = {'linestyle': '--', 'lw': 0.5}

plt.show()

plt.savefig('p')
