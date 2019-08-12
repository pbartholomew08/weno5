import math
import numpy as np
import matplotlib.pyplot as plt

import weno
weno5 = weno.weno.weno5

N = 100
L = 2 * math.pi

dx = L / (N - 1.0)
x = []
f = []
fp = []
for i in range(N):
  x.append(i * dx)
  f.append(math.sin(x[i]))
  fp.append(math.cos(x[i]))

# Test x
u = np.zeros((N, 1, 1), dtype=np.float64, order="F")
phi = np.zeros((N, 1, 1), dtype=np.float64, order="F")
gradphi = np.zeros((N, 1, 1), dtype=np.float64, order="F")
for i in range(N):
  for j in range(1):
    for k in range(1):
      u[i][j][k] = 1.0
      phi[i][j][k] = f[i]
      gradphi[i][j][k] = 0.0
weno5(gradphi, phi, u, 1, 2, 2, dx, dx, dx)
fpc = np.zeros(N)
for i in range(N):
  fpc[i] = gradphi[i][0][0]
plt.plot(x, fpc, marker="o")
plt.plot(x, fp)
plt.title("Test x-derivative (smooth)")
plt.savefig("weno-smoothx.eps", bbox_inches="tight")
plt.close()

# Test y
u = np.zeros((1, N, 1), dtype=np.float64, order="F")
phi = np.zeros((1, N, 1), dtype=np.float64, order="F")
gradphi = np.zeros((1, N, 1), dtype=np.float64, order="F")
for i in range(1):
  for j in range(N):
    for k in range(1):
      u[i][j][k] = 1.0
      phi[i][j][k] = f[j]
      gradphi[i][j][k] = 0.0
weno5(gradphi, phi, u, 2, 2, 2, dx, dx, dx)
fpc = np.zeros(N)
for i in range(N):
  fpc[i] = gradphi[0][i][0]
plt.plot(x, fpc, marker="o")
plt.plot(x, fp)
plt.title("Test y-derivative (smooth)")
plt.savefig("weno-smoothy.eps", bbox_inches="tight")
plt.close()

# Test z
u = np.zeros((1, 1, N), dtype=np.float64, order="F")
phi = np.zeros((1, 1, N), dtype=np.float64, order="F")
gradphi = np.zeros((1, 1, N), dtype=np.float64, order="F")
for i in range(1):
  for j in range(1):
    for k in range(N):
      u[i][j][k] = 1.0
      phi[i][j][k] = f[k]
      gradphi[i][j][k] = 0.0
weno5(gradphi, phi, u, 3, 2, 2, dx, dx, dx)
fpc = np.zeros(N)
for i in range(N):
  fpc[i] = gradphi[0][0][i]
plt.plot(x, fpc, marker="o")
plt.plot(x, fp)
plt.title("Test z-derivative (smooth)")
plt.savefig("weno-smoothz.eps", bbox_inches="tight")
plt.close()

# Test with discontinuity
for i in range(N/2, N):
  f[i] += 1
# Test x
u = np.zeros((N, 1, 1), dtype=np.float64, order="F")
phi = np.zeros((N, 1, 1), dtype=np.float64, order="F")
gradphi = np.zeros((N, 1, 1), dtype=np.float64, order="F")
for i in range(N):
  for j in range(1):
    for k in range(1):
      u[i][j][k] = 1.0
      phi[i][j][k] = f[i]
      gradphi[i][j][k] = 0.0

weno5(gradphi, phi, u, 1, 2, 2, dx, dx, dx)

fpc = np.zeros(N)
for i in range(N):
  fpc[i] = gradphi[i][0][0]
plt.plot(x, fpc, marker="o")
plt.title("Test x-derivative (discontinuous)")
plt.savefig("weno-discontinuousx.eps")
plt.close()
