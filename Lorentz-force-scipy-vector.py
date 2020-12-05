#Trajectory of charged particle in E & M fields, use vector equations
from pylab import *
from scipy import integrate

m = 25.0	# Mass and 
q = 5.0		# Charge of the particle
E = array([0, 0, .1]) # Electric field components Ex,Ey & Ez 
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

# Set initial values
xs[0], ys[0], zs[0] = (0., 1., 1.05)

# Step through "time", calculating the partial derivatives at the current point
# and using them to estimate the next point
for i in range(num_steps):
    x_dot, y_dot, z_dot = lorenz(xs[i], ys[i], zs[i])
    xs[i + 1] = xs[i] + (x_dot * dt)
    ys[i + 1] = ys[i] + (y_dot * dt)
    zs[i + 1] = zs[i] + (z_dot * dt)

B = np.array(xs,ys,zs)
def solver(X, t0): # X is a six element array, t0 required for the solver
	v = array([X[3], X[4], X[5]])  # make the velocity vector
	a = q * (E + cross(v,B)) / m   # F = q v x B ; a = F/m
	return [X[3], X[4], X[5], a[0], a[1], a[2]]	


tmax = 50  # calculate up to 50 seconds			
x0 = [0, 0, 0, 0, 1, 0]					# position & velocity, at t = 0
t = arange(0, tmax, 0.01)	   			# array of time coordinate
pv = integrate.odeint(solver, x0, t)	# integrate for position and velocity

savetxt('xyz.txt',pv)		# saves 6 columns to file, x,y,x, Vx, Vy, Vz

from mpl_toolkits.mplot3d import Axes3D
ax = Axes3D(figure())
ax.plot(pv[:,0], pv[:,1], pv[:,2])			# 3d plot of x, y and z
ax.set_zlabel('Z axis')
show()

