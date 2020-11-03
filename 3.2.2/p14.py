import numpy as np
import matplotlib.pyplot as plt
import math

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

data = []
with open("data-1.txt") as f:
    for line in f:
        data.append([float(x) for x in line.split()])

f1 = [0 for i in range(len(data))]
phase1 = [0 for i in range(len(data))]

for i in range(len(data)):
    f1[i] = data[i][0]
    phase1[i] = data[i][1]

data = []
with open("data-2.txt") as f:
    for line in f:
        data.append([float(x) for x in line.split()])

f2 = [0 for i in range(len(data))]
phase2 = [0 for i in range(len(data))]

for i in range(len(data)):
    f2[i] = data[i][0]
    phase2[i] = data[i][1]

plt.scatter(f2, phase2, s=5, color='navy')
plt.scatter(f1, phase1, s=5, color='orange')

gamma1 = 3.5389
gamma2 = 2.3825

t1 = np.arange(12.01, 28.49, 0.02)
q1 = 1 / 3.14 * t1 * gamma1 / (28.5 ** 2 - t1 ** 2)
for n in range(len(q1)):
    q1[n] = 1 / 3.14 * math.atan(q1[n])

t3 = np.arange(28.51, 36., 0.02)
q3 = 1 / 3.14 * t3 * gamma1 / (28.5 ** 2 - t3 ** 2)
for n in range(len(q3)):
    q3[n] = 1 + 1 / 3.14 * math.atan(q3[n])

t2 = np.arange(12., 21.49, 0.02)
q2 = 1 / 3.14 * t2 * gamma2 / (21.5 ** 2 - t2 ** 2)
for n in range(len(q2)):
    q2[n] = 1 / 3.14 * math.atan(q2[n])

t4 = np.arange(21.51, 36, 0.02)
q4 = 1 / 3.14 * t4 * gamma2 / (21.5 ** 2 - t4 ** 2)
for n in range(len(q4)):
    q4[n] = 1 + 1 / 3.14 * math.atan(q4[n])

plt.plot(t1, q1, lw=1.0, color='#36B9DA', label=r'$C_2 = 33{,}2$, нФ')
plt.plot(t3, q3, lw=1.0, color='#36B9DA')
plt.plot(t2, q2, lw=1.0, color='#B67719', label=r'$C_4 = 57{,}2$, нФ')
plt.plot(t4, q4, lw=1.0, color='#B67719')

plt.axvline(x=28.5, ymin=0., ymax=1, ls='--', color='navy', lw=1.0)

plt.axvline(x=21.5, ymin=0., ymax=1, ls='--', color='orange', lw=1.0)

plt.legend()


plt.xlabel(r'$\nu$, кГц', size=18)
plt.ylabel(r'$\Delta \varphi / \pi$', size=18)
kwargs = {'linestyle': '--', 'lw': 0.5}
plt.grid(True, **kwargs)

plt.show()

