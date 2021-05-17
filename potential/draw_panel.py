import matplotlib.pyplot as plt
import numpy as np
import lj_potential as lj

x0 = np.zeros([5 * 3])
x0 = x0 -3 * np.random.rand(x0.size).reshape(x0.shape)

if __name__ == '__main__':
    data = lj.mont_carlo(lj.lennard_jones_potential, x0, lj.first_derivative, trash=100)[0].reshape([5, 3])
    x = np.array([])
    y = np.array([])
    z = np.array([])
    print(data)
    for i in range(0, len(data)):
        x = np.append(x, data[i, 0])
        y = np.append(y, data[i, 1])
        z = np.append(z, data[i, 2])

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x, y, z, zdir='z')
    plt.show()

