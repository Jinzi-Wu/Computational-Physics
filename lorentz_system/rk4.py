import numpy as np


def runge_kutta_4th(fx, fy, fz, coor0, step, thrash=6 * 1e5):
    """
    use runge-kutta 4th method calculate the lorentz system
    :param fx: x_prime of t
    :param fy: y_prime of t
    :param fz: z_prime of t
    :param coor0: initial point
    :param step: step len
    :param thrash: total steps
    :return:
    """

    x = np.append(np.array([]), coor0[0])
    y = np.append(np.array([]), coor0[1])
    z = np.append(np.array([]), coor0[2])
    t = np.append(np.array([]), 0)
    order = 0

    while order <= thrash:
        temp_coor = coor0
        fp1 = np.array([fx(coor0[0], coor0[1]), fy(coor0[0], coor0[1], coor0[2]), fz(coor0[0], coor0[1], coor0[2])])
        coor0 = coor0 + np.array([step / 2 * fp1[0], step / 2 * fp1[1], coor0[2] + step / 2 * fp1[2]])

        fp2 = np.array([fx(coor0[0], coor0[1]), fy(coor0[0], coor0[1], coor0[2]), fz(coor0[0], coor0[1], coor0[2])])
        coor0 = coor0 + np.array([step / 2 * fp2[0], step / 2 * fp2[1], coor0[2] + step / 2 * fp2[2]])

        fp3 = np.array([fx(coor0[0], coor0[1]), fy(coor0[0], coor0[1], coor0[2]), fz(coor0[0], coor0[1], coor0[2])])
        coor0 = coor0 + np.array([step / 2 * fp3[0], step / 2 * fp3[1], coor0[2] + step / 2 * fp3[2]])

        fp4 = np.array([fx(coor0[0], coor0[1]), fy(coor0[0], coor0[1], coor0[2]), fz(coor0[0], coor0[1], coor0[2])])

        fp = fp1 + 2 * fp2 + 2 * fp3 + fp4
        coor0 = temp_coor + step / 6 * fp
        x = np.append(x, coor0[0])
        y = np.append(y, coor0[1])
        z = np.append(z, coor0[2])
        order += 1
        t = np.append(t, order*step)

    return x, y, z


def ordinary(fx, fy, fz, coor0, step, thrash=1e5):
    x = np.append(np.array([]), coor0[0])
    y = np.append(np.array([]), coor0[1])
    z = np.append(np.array([]), coor0[2])
    order = 0

    while order <= thrash:
        coor0 = coor0 + step * np.array(
            [fx(coor0[0], coor0[1]), fy(coor0[0], coor0[1], coor0[2]), fz(coor0[0], coor0[1], coor0[2])])
        x = np.append(x, coor0[0])
        y = np.append(y, coor0[1])
        z = np.append(z, coor0[2])
        order += 1
    return x, y, z
