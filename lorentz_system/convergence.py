import rk4
import numpy as np
import matplotlib.pyplot as plt

x0 = np.zeros([38 * 3])
x0 = x0 + 3.5 * np.random.rand(x0.size).reshape(x0.shape)


def fx(x, y, sigma=10):
    return sigma * (y - x)


def fy(x, y, z, rou=28):
    return x * (rou - z) - y


def fz(x, y, z, beta=8/3):
    return x * y - beta * z


coor0 = np.zeros([3])
coor0[1] += 1
m = 1
energy = []

data = rk4.ordinary(fx, fy, fz, coor0, 0.01, 1e4)



file_handle = open('coor2.txt', mode='w+')
file_handle.write('x :' + str(data[0][600 : 1000]))
file_handle.write('y :' + str(data[1][600 : 1000]))
file_handle.write('z :' + str(data[2][600 : 1000]))
file_handle.close()


x = data[0]
y = data[1]
z = data[2]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(x, y, z, 'b.')
ax.legend('sigma=10, rou=28, beta=8/3')
plt.show()

print(data[0][600 : 1000])
