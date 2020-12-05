# Charged particle trajectory under electric and magnetic fields
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import seaborn as sbn
m = 25.0  # Mass and
q = 5.0		# Charge of the particle
Ex = 0.0  # Electric Field vector
Ey = 0.0
Ez = 0.1
Bx = 0.0  # Magnetic field vector
By = 0.0
Bz = 5.0


def lorenz(x, y, z, s=10, r=28, b=2.667):
    '''
    Given:
       x, y, z: a point of interest in three dimensional space
       s, r, b: parameters defining the lorenz attractor
    Returns:
       x_dot, y_dot, z_dot: values of the lorenz attractor's partial
           derivatives at the point x, y, z
    '''
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot


dt = 0.01
num_steps = 10000

# Need one more for the initial values
xs = np.empty(num_steps + 1)
ys = np.empty(num_steps + 1)
zs = np.empty(num_steps + 1)
ax = np.empty(num_steps + 1)
ay = np.empty(num_steps + 1)
az = np.empty(num_steps + 1)

# Set initial values
xs[0], ys[0], zs[0] = (0., 1., 1.05)
vx = 0.1
vy = 0.0
vz = 0.0
# Step through "time", calculating the partial derivatives at the current point
# and using them to estimate the next point
for i in range(num_steps):
    x_dot, y_dot, z_dot = lorenz(xs[i], ys[i], zs[i])
    xs[i + 1] = (xs[i] + (x_dot * dt))
    ys[i + 1] = (ys[i] + (y_dot * dt))
    zs[i + 1] = (zs[i] + (z_dot * dt))
    ax[i] = q * (Ex + (vy * zs[i]) - (vz * ys[i])) / m  # Lorentz force / mass
    ay[i] = q * (Ey - (vx * zs[i]) + (vz * xs[i])) / m
    az[i] = q * (Ez + (vx * ys[i]) - (vy * xs[i])) / m

ax = np.array(ax)
ay = np.array(ay)
az = np.array(az)

import matplotlib.pyplot as plt
plt.scatter(ax,ay)
plt.show()
plt.scatter(ay,az)
plt.show()
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(ax, ay, az)
plt.show()