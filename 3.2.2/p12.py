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
    #1 approximation : y = kx + b
    if (d == 1):
        for n in range(len(X)):
            xyav = xyav + X[n]*Y[n] / m
            y2av = y2av + Y[n]*Y[n] / m
            x2av = x2av + X[n]*X[n] / m
            xav = xav + X[n]/m
            yav = yav + Y[n]/m
        k = ((xyav) - xav * yav)/(x2av - xav*xav)
        b = yav - k * xav
        d_k = 1/math.sqrt(m) * math.sqrt(abs((y2av - yav * yav)/(x2av - xav * xav) - k**2))
        d_b = d_k * math.sqrt(x2av - xav*xav)
        final_data = [k, d_k, b, d_b]
    else:
    #2 approximation : y = kx
        if (d == 2):
            for n in range(len(X)):
                xyav = xyav + X[n]*Y[n] / m
                y2av = y2av + Y[n]*Y[n] / m
                x2av = x2av + X[n]*X[n] / m
            k = (xyav)/(x2av)
            d_k = 1/math.sqrt(m) * math.sqrt (abs((y2av)/(x2av) - k**2))
            final_data = [k, d_k]
    return final_data

X1 = [26.2, 27.2, 28., 29.3, 30., 31.4, 32.1, 33, 34.3, 35.2, 36.3, 25, 23.7, 22.8, 21.4, 20]
Y1 = [3.1, 6.0, 6.4, 3.14, 2.2, 1.38, 1.14, 0.93, 0.74, 0.63, 0.54, 1.81, 1.31, 1.1, 0.89, 0.76]

gamma1 = 3.5389

X2 = [20.2, 21., 22.2, 23.7, 24.8, 25.6, 26.8, 27.3, 28, 19.5, 18.1, 17.2, 16.4, 15.8, 14.3, 13.7]
Y2 = [3.5, 6.56, 3.39, 1.47, 1.01, 0.82, 0.64, 0.59, 0.51, 2.23, 1.34, 1.05, 0.91, 0.82, 0.68, 0.63]

gamma2 = 2.3825

X_1 = [0 for i in range(len(X1))]

for i in range(len(X1)):
        X_1[i] = 1 / math.sqrt((3.14**2)*(28.5**2 - X1[i]**2)**2 + (gamma1**2)*(X1[i]**2))

X_2 = [0 for i in range(len(X2))]

for i in range(len(X2)):
        X_2[i] = 1 / math.sqrt((3.14**2)*(21.5**2 - X2[i]**2)**2 + (gamma2**2)*(X2[i]**2))



data = approx(X_1, Y1, 2)

k1 = data[0]
d_k1 = data[1]

data = approx(X_2, Y2, 2)

k2 = data[0]
d_k2 = data[1]

plt.scatter(X1, Y1, s=5, color='navy')
plt.scatter(X2, Y2, s=5, color='orange')


t1 = np.arange(12, 36, 0.1)
q1 = 9.86 * (28.5**2 - t1 * t1)**2 + gamma1**2 *(t1 * t1)
s1 = k1*k1 / q1
for n in range(len(s1)):
        s1[n] = math.sqrt(s1[n])

t2 = np.arange(12, 36, 0.1)
q2 = 9.86 * (21.5**2 - t2 * t2)**2 + gamma2**2 *(t2 * t2)
s2 = k2*k2 / q2
for n in range(len(s2)):
        s2[n] = math.sqrt(s2[n])       

plt.plot(t1, s1, lw=1.0, color='#36B9DA', label=r'$C_2 = 33{,}2$, нФ')
plt.plot(t2, s2, lw=1.0, color='#B67719', label=r'$C_4 = 57{,}2$, нФ')


plt.axvline(x=28.5, ymin=0., ymax=1.8, ls='--', color='navy', lw=1.0)
plt.axvline(x=21.5, ymin=0., ymax=1.8, ls='--', color='orange', lw=1.0)
plt.legend(fontsize='x-large')


plt.xlabel(r'$\nu$, кГц', size=18)
plt.ylabel(r'$U_C$, В', size=18)
kwargs={'linestyle':'--', 'lw':0.5}
plt.grid(True, **kwargs)
plt.legend()
plt.show()

plt.savefig(r'p12.png')

