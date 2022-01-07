from scipy.optimize import fsolve
from scipy.integrate import quad
import numpy as np


def g(y, z):
    return np.cos(y - z)


def f(x):
    eq = -1 + quad(g, np.pi / 2, np.pi, args=(x,))[0]
    return eq

#res, err = quad(g, np.pi / 2, np.pi, args=(0 ,))
res = fsolve(f, 1.5)

print(res)
print(np.pi/2)
