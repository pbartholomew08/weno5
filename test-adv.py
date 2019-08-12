from math import sin, pi, log, log10, sqrt, exp
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import weno
weno5 = weno.weno.weno5

def init_pulse(x, xl, xr):

    phi = []
    for i in range(len(x)):
        if (x[i] >= xl) and (x[i] <= xr):
            xp = x[i] - xl

            # phi.append(sin((xp/(xr - xl)) * pi))
            phi.append(1)
        else:
            phi.append(0)

    return phi
def init_jiang(x):

    phi = []
    n = len(x)

    a = 0.5
    z = -0.7
    d = 0.005
    alpha = 10.0
    beta = log10(2.0) / (36 * d**2)

    for i in range(n):
        if (-0.8 <= x[i]) and (x[i] <= -0.6):
            phi.append(g(x[i], beta, z - d) + g(x[i], beta, z + d) + 4 * g(x[i], beta, z))
            phi[-1] /= 6.0
        elif (-0.4 <= x[i]) and (x[i] <= -0.2):
            phi.append(1)
        elif (0 <= x[i]) and (x[i] <= 0.2):
            phi.append(1 - abs(10 * (x[i] - 0.1)))
        elif (0.4 <= x[i]) and (x[i] <= 0.6):
            phi.append(f(x[i], alpha, a - d) + f(x[i], alpha, a + d) + 4 * f(x[i], alpha, a))
            phi[-1] /= 6.0
        else:
            phi.append(0)

    return phi

def g(x, b, z):
    return exp(-b * (x - z)**2)

def f(x, alpha, a):
    return sqrt(max(1 - (alpha**2) * (x - a)**2, 0))
def calc_rhs(t, y, f_args):
    u = f_args[0]  # The velocity field
    dx = f_args[1] # The grid spacing
    n = len(y)

    y3d = np.array(y).reshape((n, 1, 1), order="F")
    u3d = u * np.ones(n).reshape((n, 1, 1), order="F")
    dydx = np.zeros((n, 1, 1), order="F")

    weno5(dydx, y3d, u3d, 1, 0, 0, dx, dx, dx)

    return -u*dydx.reshape(n)

class rk3():

    def __init__(self, f, t = 0):

        self.f = f
        self.t = 0

    def set_initial_value(self, y0):

        self.y = y0

    def set_f_params(self, f_args):
        self.f_args = f_args

    def successful(self):
        return True

    def integrate(self, tnext):

        dt = tnext - self.t

        # Stage 1
        f0 = self.f(self.t, self.y, self.f_args)
        y1 = self.y + dt * f0

        # Stage 2
        f1 = self.f(self.t, y1, self.f_args)
        y2 = self.y + (dt / 4.0) * (f0 + f1)

        # Stage 3
        f2 = self.f(self.t, y2, self.f_args)
        self.y += (dt / 6.0) * (f0 + 4 * f2 + f1)

        self.t += dt

L=2.0
U=1.0
N=200
CFL = 0.2
T=10 #10.0 * (L / U)

dx=L/float(N)
x = []
for i in range(N):
    x.append(i * dx - 1)
xl = -0.2
xr = 0.2

dt = CFL * dx / U

# r = ode(calc_rhs).set_integrator("dopri5", atol=1.0e-16, rtol=1.0e-8)
r = rk3(calc_rhs)
# r.set_initial_value(init_pulse(x, xl, xr))
r.set_initial_value(init_jiang(x))
r.set_f_params((U, dx))

passed_eight = False
while r.successful() and r.t < T:
    if r.t == 0:
        plt.plot(x, r.y, color="black")
    elif (r.t >= 8) and (not passed_eight):
        plt.plot(x, r.y, ls="", marker="o", color="blue")
        passed_eight = True
    print r.t, min(r.y), max(r.y)
    r.integrate(r.t+dt)

plt.plot(x, r.y, ls="", marker="o", color="red")
plt.savefig("adv_test.eps", bbox_inches="tight")
