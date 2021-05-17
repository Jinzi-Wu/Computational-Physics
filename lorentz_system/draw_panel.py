import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

import rk4


def fx(x, y, sigma=10):
    return sigma * (y - x)


def fy(x, y, z, rou=50):
    return x * (rou - z) - y


def fz(x, y, z, beta=8/3):
    return x * y - beta * z


coor0 = np.zeros([3])
coor0[1] += 1
m = 1
energy = []

data = rk4.runge_kutta_4th(fx, fy, fz, coor0, 0.01, 1e4)
x = data[0]
y = data[1]
z = data[2]
t = np.linspace(0, 100, 10000+2)
for xs, ys, zs in zip(x, y, z):
   energy = np.append(energy, m * (xs**2 + ys**2 + zs**2)/2)


fig = plt.figure()
fig1, ax1 = plt.subplots(1, 1)
fig2, ax2 = plt.subplots(1, 1)
ax = fig.add_subplot(projection='3d')



ax.plot(x, y, z, '')
ax.legend('sigma=10, rou=28, beta=8/3')

plt.title('lorentz system ')
plt.savefig('lorentz system rou=1.png')

ax1.plot(x, z, 'g-')
plt.xlabel('x')
plt.ylabel('z')
plt.savefig('lorentz system on x-z.png')

ax2.plot(t, energy, 'y-')
plt.xlabel('t')
plt.ylabel('energy')
plt.savefig('lorentz system energy variation.png')


plt.show()
