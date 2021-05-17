import math


def g(x):
    return (x + 1) ** (1 / 5)


def potential(x, a=1, k=1):
    return (k/x)*math.exp(-a*x)


def f(x, b=1, energy=1):
    return 1 - (b / x) ** 2 - potential(x) / energy


def iteration(x0, thrash):
    s = 1
    while abs((x1 := g(x0)) - x0) > thrash:
        x0 = x1
        s += 1
    return x1, s


def bisection(x1, x2, thrash, b):
    while abs(x1 - x2) >= thrash:
        y1 = f(x1, b)
        y2 = f(x2, b)
        if y1 * y2 < 0:
            xm = (x1 + x2) / 2
            ym = f(xm, b)
            if y1 * ym < 0:
                x2 = xm
            elif y2 * ym < 0:
                x1 = xm
            elif ym == 0:
                return xm
    return (x1 + x2) / 2


def secant_method(x1, x2, thrash):
    while abs(x1 - x2) > thrash:
        x0 = x2 - f(x2) * (x2 - x1) / (f(x2) - f(x1))
        if f(x0) == 0:
            return x0
        elif f(x1) == 0:
            return x1
        elif f(x2) == 0:
            return x2
        if f(x1) * f(x0) < 0:
            x2 = x0
            print("x2 = " + str(x2))
        elif f(x2) * f(x0) < 0:
            x1 = x0
            print("x1 = " + str(x1))
        # print(str(x1) + "  " + str(x2))
    return x2 - f(x2) * (x2 - x1) / (f(x2) - f(x1))


if __name__ == '__main__':
    print(bisection(0.01, 100, 1e-4, 10))
