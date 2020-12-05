# polar to cartesian
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 
Bx = 200
By = 20
Bz = 20
vx = 1
vy = 1
vz = 1
X = 1.0
Y = 1.0 
Z = 3.0 
K = 200
x1 = []
y1 = []
z1 = []
dt = 0.0001
nsteps = 20000
for i in range(nsteps):
    r = np.sqrt(X*X+Y*Y+Z*Z)
    theta = np.arccos(Z/r)
    phi = np.arctan(Y/X)

    # magnetic field/mass
    Bx = K*(1/(r**(3)))*(2*np.cos(theta)*np.sin(theta)*np.cos(phi)+np.sin(theta)*np.cos(theta)*np.cos(phi))
    By = K*(1/(r**(3)))*(2*np.cos(theta)*np.sin(theta)*np.sin(phi)+np.sin(theta)*np.cos(theta)*np.sin(phi))
    Bz = K*(1/(r**(3)))*(2*np.cos(theta)*np.cos(theta)-np.sin(theta)*np.sin(theta))

    # acceleration components
    ax = (1.6*10**(-19))*(vy*Bz - vz*By)
    ay = (1.6*10**(-19))*(vz*Bx - vx*Bz)
    az = (1.6*10**(-19))*(vx*By - vy*Bx)

    # velocity components
    vx = vx + ax*dt
    vy = vy + ay*dt
    vz = vz + az*dt

    # position components
    X = X + vx*dt + 0.5*ax*dt*dt
    Y = Y + vy*dt + 0.5*ay*dt*dt
    Z = Z + vz*dt + 0.5*az*dt*dt

    # add position values to position vectors
    x1.append(X)
    y1.append(Y)
    z1.append(Z)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x1, y1, z1)
plt.show()