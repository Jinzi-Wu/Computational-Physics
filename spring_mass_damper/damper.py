import numpy as np


def inverse(A):
    """
    to calculate the inverse of a 2×2 matrix
    :param A: original matrix
    :return: the inverse of A
    """
    a, b, c, d = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
    a_inverse = np.zeros([2, 2])
    if a == 0:
        if c == 0:
            raise ValueError('invalid update matrix')
        elif c != 0:
            a_inverse[0, 0], a_inverse[0, 1], a_inverse[1, 0], a_inverse[1, 1] = d / (a * d - b * c), b / (
                        b * c - a * d), c / (b * c - a * d), a / (a * d - b * c)
    elif a != 0:
        a_inverse[0, 0], a_inverse[0, 1], a_inverse[1, 0], a_inverse[1, 1] = d / (
                -b * c + a * d), - b / (d * a - b * c), c / (b * c - a * d), a / (a * d - b * c)
    print(a_inverse)
    return a_inverse


def forward_euler(y0, A, step, thrash=1e5):
    """

    :param step: step len
    :param y0: the initial point of the system
    :param A: update matrix
    :param thrash: update steps
    :return:
    """
    sp = 0
    t = 0
    vx_x = np.append(np.array([]), y0[0, 0])
    vx_v = np.append(np.array([]), y0[1, 0])
    xt_t = np.append(np.array([]), t)
    xt_x = np.append(np.array([]), y0[0, 0])
    a_inverse = inverse(A)
    while sp <= thrash:
        y_prime = a_inverse @ y0
        y0 = y0 + y_prime * step
        vx_x = np.append(vx_x, y0[0, 0])
        vx_v = np.append(vx_v, y0[1, 0])
        xt_t = np.append(xt_t, t)
        xt_x = np.append(xt_x, y0[0, 0])
        sp += 1
        t += step
    print(vx_x)
    print(vx_v)
    print(xt_x)
    print(xt_t)
    return vx_x, vx_v, xt_t, xt_x


def backward_euler(y0, A, step, thrash=1e5):
    """
    use euler backward method
    :param y0: initial point of the system
    :param A: prime relation
    :param step: step len
    :param thrash: total steps
    :return:
    """
    sp = 0
    t = 0
    vx_x = np.append(np.array([]), y0[0, 0])
    vx_v = np.append(np.array([]), y0[1, 0])
    xt_t = np.append(np.array([]), t)
    xt_x = np.append(np.array([]), y0[0, 0])
    a_inverse = inverse(np.diag([1, 1])-step * inverse(A))
    while sp <= thrash:
        y_prime = a_inverse @ y0
        y0 = y0 + y_prime * step
        vx_x = np.append(vx_x, y0[0, 0])
        vx_v = np.append(vx_v, y0[1, 0])
        xt_t = np.append(xt_t, t)
        xt_x = np.append(xt_x, y0[0, 0])
        sp += 1
        t += step
    print(vx_x)
    print(vx_v)
    print(xt_x)
    print(xt_t)
    return vx_x, vx_v, xt_t, xt_x

