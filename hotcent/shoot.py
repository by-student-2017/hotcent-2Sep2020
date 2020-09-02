import numpy as np


def shoot(u, dx, c2, c1, c0, N):
    """ Integrate diff equation

           2
         d u      du
         --- c  + -- c  + u c  = 0
           2  2   dx  1      0
         dx

    in equispaced grid (spacing dx) using simple finite difference formulas

    u'(i) = (u(i+1) - u(i-1)) / (2*dx) and
    u''(i) = (u(i+1) - 2*u(i) + u(i-1)) / dx**2

    u[0:2] *has already been set* according to boundary conditions.

    return u, number of nodes, the discontinuity of derivative at
    classical turning point (ctp), and ctp
    c0(r) is negative with large r, and turns positive at ctp.
    """
    fp = c2 / dx ** 2 + 0.5 * c1 / dx
    fm = c2 / dx ** 2 - 0.5 * c1 / dx
    f0 = c0 - 2 * c2 / dx ** 2
 
    # backward integration down to classical turning point ctp
    # (or one point beyond to get derivative)
    # If no ctp, integrate half-way
    u[-1] = 1.0
    u[-2] = u[-1] * f0[-1] / fm[-1]
    all_negative = np.all(c0 < 0)
    for i in range(N - 2 , 0, -1):
        u[i - 1] = (-fp[i] * u[i + 1] - f0[i] * u[i]) / fm[i]
        if abs(u[i - 1]) > 1e10:
            u[i - 1:] *= 1e-10  # numerical stability
        if c0[i] > 0:
            ctp = i
            break
        if all_negative and i == N // 2:
            ctp = N // 2
            break

    utp = u[ctp]
    utp1 = u[ctp + 1]
    dright = (u[ctp + 1] - u[ctp - 1]) / (2 * dx)

    for i in range(1, ctp + 1):
        u[i + 1] = (-f0[i] * u[i] - fm[i] * u[i - 1]) / fp[i]

    dleft = (u[ctp + 1] - u[ctp - 1]) / (2 * dx)
    scale = utp / u[ctp]
    u[:ctp + 1] *= scale
    u[ctp + 1] = utp1  # above overrode
    dleft *= scale
    u = u * np.sign(u[1])

    nodes = np.sum((u[0:ctp - 1] * u[1:ctp]) < 0)
    A = (dright - dleft) * utp
    return u, nodes, A, ctp
