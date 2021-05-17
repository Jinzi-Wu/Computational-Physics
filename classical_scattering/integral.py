import math

import numpy as np

import find_root as fd


def func1(x, b):
    return 1 / (x ** 2 * (1 - (b / x) ** 2) ** 0.5)


def func2(x, b):
    a = 1 / (x ** 2 * math.sqrt(1 - (b / x) ** 2 - (1 / x) * math.exp(-x)))
    return a


def integral(lower, h, b):
    print('in integral lower = ' + str(lower))
    x01 = b + h/5
    x02 = lower + h
    s1 = 0
    s2 = 0
    while x01 < 1000 and x02 < 1000:
        s1 += h / 3 * func1(x01 , b) + 4 * h / 3 * func1(x01 + h, b) + h / 3 * func1(x01 + 2*h, b)
        x01 += h * 2
        s2 += h / 3 * func2(x02 , b) + 4 * h / 3 * func2(x02 + h, b) + h / 3 * func2(x02 + 2*h, b)
        x02 += h * 2
    return s1, s2


if __name__ == '__main__':
    def f(b):
        lower = fd.bisection(1e-2, 10, 1e-4, b)
        itl = integral(lower, 1e-3, b)
        a = 2 * b * (itl[0] - itl[1]) - math.pi / 2
        print(lower)
        print(a)
        return a


    def bisection(x1, x2, thrash):
        while abs(x1 - x2) >= thrash:
            xm = (x1 + x2) / 2
            print('x1 = ' + str(x1) + ' x2 = ' + str(x2) + ' xm = ' + str(xm))
            y1 = f(x1)
            y2 = f(x2)
            ym = f(xm)
            if y1 * y2 < 0:
                if y1 * ym < 0:
                    x2 = xm
                    print('x2 = xm')
                elif y2 * ym < 0:
                    x1 = xm
                    print('x1 = xm')
                elif ym == 0:
                    return xm
            elif y1 == 0:
                return x1
            elif y2 == 0:
                return x2
            else:
                print(y1 * y2)
                raise ValueError('No root in this interval, please choose another beginning interval!')
        return (x1 + x2) / 2


    print(bisection(0.1, 1, 1e-4))
    print('left = ' + str(f(0.1)))
    print('right = ' + str(f(1)))
    # i = np.logspace(-10, 2, 10, base=2)
    # y = np.array([f(x)
    #               for x in i])
    # print(i, y)
