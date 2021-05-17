import math
import matplotlib.pyplot as plt
import numpy as np
import numpy as np


def func1(x, b):
    return 1 / (x ** 2 * (1 - (b / x) ** 2) ** 0.5)


def func2(x, b, ep):
    a = 1 / (x ** 2 * (1 - (b / x) ** 2 - potential(x, ep)) ** 0.5)
    return a


def reverse(b, ag, ep):
    lower = bisection(0.001, 1000, 1e-4, b, ep)
    print('lower = ' + str(lower))
    a = 2 * b * integral(lower, 1e-2, b, ep) - ag
    print('a = ' + str(a))
    return a


def integral(lower, h, b, ep):
    x01 = b + h
    x02 = lower + h
    s1 = 0
    s2 = 0
    print('b in integral ' + str(b) + ' b/x = ' + str(b / x01))
    while x01 < 1000 and x02 < 1000:
        s1 += h / 3 * func1(x01, b) + 4 * h / 3 * func1(x01 + h, b) + h / 3 * func1(x01 + 2 * h, b)
        x01 += h * 2
        s2 += h / 3 * func2(x02, b, ep) + 4 * h / 3 * func2(x02 + h, b, ep) + h / 3 * func2(x02 + 2 * h, b, ep)
        x02 += h * 2
    print(s1 - s2)
    return s1 - s2


def potential(r, a):
    r6 = (1 / r) ** 6

    return 4 * a * r6 * (r6 - 1)


def f(x, ep, b=1, energy=1):
    return 1 - (b / x) ** 2 - potential(x, ep) / energy


def bisection(x1, x2, thrash, b, ep):
    while abs(x1 - x2) >= thrash:
        y1 = f(x1, ep, b)
        y2 = f(x2, ep, b)
        if y1 * y2 < 0:
            xm = (x1 + x2) / 2
            ym = f(xm, ep, b)
            if y1 * ym < 0:
                x2 = xm
            elif y2 * ym < 0:
                x1 = xm
            elif ym == 0:
                return xm
    return (x1 + x2) / 2


def bisection2(x1, x2, thrash, ag, ep):
    print('ag = ' + str(ag))
    while abs(x1 - x2) >= thrash:
        y1 = reverse(x1, ag, ep)
        y2 = reverse(x2, ag, ep)
        if y1 * y2 < 0:
            xm = (x1 + x2) / 2
            ym = reverse(xm, ag, ep)
            if y1 * ym < 0:
                x2 = xm
            elif y2 * ym < 0:
                x1 = xm
            elif ym == 0:
                return xm
    return (x1 + x2) / 2


def differential_cross_section(b, ep):
    print('b = ' + str(b))
    lower = bisection(0.001, 1000, 1e-5, b, ep)
    theta1 = 2 * b * integral(lower, 0.01, b - 1e-8, ep)
    theta2 = 2 * b * integral(lower, 0.01, b + 1e-8, ep)
    a = (b / math.sin(integral(lower, 0.01, b, ep))) * abs(2 * 1e-8 / (theta2 - theta1))
    print('point value = ' + str(a))
    return a


if __name__ == '__main__':
    angle1 = np.logspace(-2, -1, 10, base=10)
    angle = np.logspace(-3, 0.82, 10, base=2)
    sp_angle = np.logspace(-2, 0.12, 10, base=2)
    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel('θ')
    ax.set_ylabel('differential cross section')
    ax.set_title('Differential Cross Section Related to θ')
    z1 = np.array([bisection2(1, 99, 1e-4, ag, 0.1)
                   for ag in angle1[:5]])
    z1 = np.append(z1, [bisection2(0.05, 10, 1e-4, ag, 1) for ag in angle[5:]])
    y1 = np.array([differential_cross_section(z, 0.1)
                   for z in z1])
    # z2 = np.array([bisection2(0.1, 99, 1e-4, ag, 1)
    #                for ag in angle1[:5]])
    # z2 = np.append(z2, [bisection2(0.1, 10, 1e-4, ag, 1)
    #                     for ag in angle[5:]])
    # y2 = np.array([differential_cross_section(z, 1)
    #                for z in z1])
    # z3 = np.array([bisection2(0.1, 99, 1e-4, ag, 10)
    #                for ag in sp_angle])
    # y3 = np.array([differential_cross_section(z, 10)
    #                for z in z3])
    z4 = np.array([bisection2(0.1, 99, 1e-4, ag, 0.01)
                   for ag in angle1[:5]])
    z4 = np.append(z4, [bisection2(0.05, 10, 1e-4, ag, 0.01)
                        for ag in angle[5:]])
    y4 = np.array([differential_cross_section(z, 0.01)
                   for z in z4])
    # y3 = np.array([1 / 16 * (1 / math.sin(ag / 2)) ** 4
    #                for ag in angle])
    # ax.plot(angle, y2, 'r.')
    # ax.plot(angle, y3)
    # plt.legend(['numerical result', 'analytical result'])
    # y2 = np.append(y2, [1 / 16 * (1 / math.sin(ag / 2)) ** 4
    #                for ag in angle[5:]])
    ax.plot(angle, y4, 'r.--')
    ax.plot(angle, y1, 'y.--')
    plt.legend(['epsilon = 0.01', 'epsilon = 0.1'])
    plt.savefig('./lj_potential_scattering_1.jpg')
    # fig, ax = plt.subplots(1, 1)
    # ax.set_xlabel('θ')
    # ax.set_ylabel('differential cross section')
    # ax.set_title('Differential Cross Section Related to θ')
    # ax.plot(angle, y2, 'g.--')
    # ax.plot(sp_angle, y3, 'b.--')
    # # print(z1)
    # # ax.plot(angle, y2)
    # plt.legend(['epsilon = 1', 'epsilon = 10'])
    # plt.savefig('./lj_potential_scattering.jpg')
    plt.show()
    # i = np.logspace(-10, 3, 100, base=2)
    # y = np.array([reverse(ia,2, 0.1)
    #               for ia in i])
    # ax.plot(i, y)
    # plt.show()
    print('left = ' + str(reverse(0.05, 1.771535038204721, 1)))
    print('right = ' + str(reverse(10, 1.771535038204721, 1)))
    # print(bisection2(0, 10, 1e-4, 2**1.65))
