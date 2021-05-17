import math
import numpy as np
import bfgs
import matplotlib.pyplot as plt
from time import time


def norm(euler_cor):
    # print('the shape of imput coordinate in norm is' + str(np.shape(euler_cor)))
    return (euler_cor[0] ** 2 + euler_cor[1] ** 2 + euler_cor[2] ** 2) ** (1 / 2)


def first_derivative(euler_coordinate, sigma=1.0, epsilon=1.0):
    """
       use five point methods to calculate the gradiant of func
       :param epsilon:
       :param sigma:
       :param euler_coordinate: array([n,1]) which is the rational coordinate of particles
       :return: the first derivative of the function array([n,1])
       """
    euler_coordinate = euler_coordinate.reshape([int(euler_coordinate.size / 3), 3])
    f_prime = np.array([])
    j = 0
    for i in range(0, len(euler_coordinate)):
        tanpor_prime = np.zeros([3])
        while j < len(euler_coordinate):
            if i != j:
                r = euler_coordinate[i] - euler_coordinate[j]
                r_norm = norm(euler_coordinate[i] - euler_coordinate[j])
                r7 = (sigma / r_norm) ** 7
                r6 = (sigma / r_norm) ** 6
                para = 24 * epsilon * r7 * (-2 * r6 + 1) * (1 / r_norm)
                tanpor_prime += r * para
            j += 1
        f_prime = np.append(f_prime, tanpor_prime)
    return f_prime


def lennard_jones_potential(euler_coor, sigma=1.0, epsilon=1.0):
    """

    Calculates the Lennard-Jones potential for particles with diameter sigma
    at a separation r with a well-depth epsilon.
    #
    # >>> lennard_jones_potential(1.0, 1.0, 1.0)
    # 0.0
    #
    # >>> lennard_jones_potential(2**(1/6), 1.0, 1.0)
    # -1.0
    #
    # >>> lennard_jones_potential(0.0, 1.0, 1.0)
    # Traceback (most recent call last):
    # ZeroDivisionError: float division by zero
    #
    # >>> lennard_jones_potential(-1.0, 1.0, 1.0)
    # Traceback (most recent call last):
    # ValueError: distance between particles is negative
    #
    # >>> lennard_jones_potential(1.0, -1.0, 1.0)
    # Traceback (most recent call last):
    ValueError: particle diameter is not strictly positive

    :param euler_coor: the euler coordinates of the particles
    :type euler_coor: array
    :param sigma: the diameter of a particle
    :type sigma: float
    :param epsilon: the well depth of the potential
    :type epsilon: float

    :return: the Lennard-Jones energy of the particle pair
    :rtype: float
    """
    euler_coor = euler_coor.reshape([int(len(euler_coor) / 3), 3])
    potential = 0
    for i in range(0, len(euler_coor)):
        for j in range(0, len(euler_coor)):
            if j != i:
                r = norm(euler_coor[i] - euler_coor[j])
                r6 = (sigma / r) ** 6
                potential += 4 * epsilon * r6 * (r6 - 1)
    return potential / 2


def mont_carlo(potential_func, x0, df, trash=1000):
    print('in mont carlo')
    coordinate = bfgs.bfgs(potential_func, x0, df)
    potential = potential_func(coordinate)
    order = 0
    recurent = 0
    t0 = time()
    while order < trash * int(math.log(int(len(x0)/3))):
        temp_coor = bfgs.bfgs(potential_func, x0 + 3.5 * np.random.rand(x0.size).reshape(x0.shape), df)
        temp_energy = potential_func(temp_coor)
        t = time() - t0
        if temp_energy < potential:
            potential, x0 = temp_energy, temp_coor
            order += 1
            print(potential)
            t0 = time()
        elif temp_energy < 0 and recurent < 100:
            recurent += 1
        elif t > 6000 * int(math.log(int(len(x0)/3))):
            return x0, potential_func(x0)
    return x0, potential_func(x0)


def mont_carlo_b(potential_func, x0, df, trash=20000):
    print('in mont carlo b')
    count = 0
    coordinate = bfgs.bfgs(potential_func, x0, df)
    potential = potential_func(coordinate)
    order = 0
    recurent = 0
    while order < trash * int(math.log(int(len(x0) / 3))):
        temp_coor = bfgs.bfgs(potential_func, x0 + 3.5 * np.random.rand(x0.size).reshape(x0.shape), df)
        temp_energy = potential_func(temp_coor)
        order += 1
        t = time() - t0
        if temp_energy < potential:
            if abs(potential - temp_energy) > 1e-5:
                potential, x0 = temp_energy, temp_coor
                count += 1
                print(potential)
                t0 = time()
        elif temp_energy < 0 and recurent < 100:
            recurent += 1
        elif t > 600 * int(math.log(int(len(x0)/3))):
            return count
    return count


if __name__ == '__main__':
    order = np.array([])
    energy = np.array([])
    counts = np.array([])

    for i in range(3, 5):
        order = np.append(order, i)
        x = np.ones([i, 3])
        y = - x + np.random.rand(x.size).reshape(x.shape)
        print('input is ' + str(y))
        result = mont_carlo(lennard_jones_potential, y.reshape([3 * len(y), ]), first_derivative)
        energy = np.append(energy, result[1])
        print('the lowest energy state is ' + '\n' +
              str(result[0].reshape((int(len(result[0]) / 3), 3))) + '\n' +
              'the lowest energy is ' + str(result[1]))

    order = np.append(order, 5)
    x = np.loadtxt('lj5.txt')
    y = x + 0.01 * np.random.rand(x.size).reshape(x.shape)
    print('input is ' + str(y))
    result = mont_carlo(lennard_jones_potential, y.reshape([3 * len(y), ]), first_derivative)
    energy = np.append(energy, result[1])
    print('the lowest energy state is ' + '\n' +
          str(result[0].reshape((int(len(result[0]) / 3), 3))) + '\n' +
          'the lowest energy is ' + str(result[1]))

    for i in range(6, 13):
        order = np.append(order, i)
        x = np.ones([i, 3])
        y = - x + 2 * np.random.rand(x.size).reshape(x.shape)
        print('input is ' + str(y))
        result = mont_carlo(lennard_jones_potential, y.reshape([3 * len(y), ]), first_derivative)
        energy = np.append(energy, result[1])
        print('the lowest energy state is ' + '\n' +
              str(result[0].reshape((int(len(result[0]) / 3), 3))) + '\n' +
              'the lowest energy is ' + str(result[1]))

    order = np.append(order, 13)
    x = np.loadtxt('_lj13.txt')
    y = x + 0.01 * np.random.rand(x.size).reshape(x.shape)
    print('input is ' + str(y))
    result = mont_carlo(lennard_jones_potential, y.reshape([3 * len(y), ]), first_derivative)
    energy = np.append(energy, result[1])
    print('the lowest energy state is ' + '\n' +
          str(result[0].reshape((int(len(result[0]) / 3), 3))) + '\n' +
          'the lowest energy is ' + str(result[1]))

    # for i in range(3, 14):
    #     x = -0.2 * np.ones([i * 3, ])
    #     y = x + 0.01 * np.random.rand(x.size).reshape(x.shape)
    #
    #     order = np.append(order, i)
    #     counts = np.append(counts, mont_carlo_b(lennard_jones_potential, y, first_derivative))
    #     print(counts)


    # plt.subplots(1, 1)
    # plt.plot(order, energy, '.-')
    # plt.legend(['the minimal energy of particles related to particle numbers'])
    # plt.ylabel('minimal energy')
    # plt.xlabel('numbers of particles')
    # plt.savefig('./MC3.png')

    plt.subplots(1, 1)
    plt.plot(order, counts, '.-')
    plt.legend(['the numbers of minima related to particle numbers'])
    plt.ylabel('numbers of minima')
    plt.xlabel('numbers of particles')
    plt.savefig('./MC4.png')
    plt.show()
