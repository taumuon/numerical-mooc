import numpy
import math

# model parameters:
g = 9.8      # gravity in m s^{-2}
v_t = 4.9  # trim velocity in m s^{-1}   
C_D = 1./5.  # drag coefficient --- or D/L if C_L=1
C_L = 1.0    # for convenience, use C_L = 1

### set initial conditions ###
v0 = v_t     # start at the trim velocity (or add a delta)
theta0 = 0.25 # initial angle of trajectory
x0 = 0.0     # horizotal position is arbitrary
y0 = 2.0  # initial altitude

T = 20.0  

def f(u):
    """Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : array of float
        array containing the solution at time n.
        
    Returns
    -------
    dudt : array of float
        array containing the RHS given u.
    """
    
    v = u[0]
    theta = u[1]
    # x = u[2]
    # y = u[3]
    return numpy.array([-g*math.sin(theta) - C_D/C_L*g/v_t**2*v**2,
                      -g*math.cos(theta)/v + g/v_t**2*v,
                      v*math.cos(theta),
                      v*math.sin(theta)])
                      
def euler_step(u, f, dt):
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equation.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    
    return u + dt * f(u)

beginTheta = 0.
endTheta = 0.5
numThetas = 10
thetas = numpy.linspace(beginTheta, endTheta, numThetas) 

x0 = 0.
y0 = 2.0

beginV0 = 1
endV0 = 1000
numVelocities = 100
velocities = numpy.linspace(beginV0, endV0, numVelocities) 

dt = 0.1

# shape is rows, columns
distances = numpy.zeros(shape=(numThetas, numVelocities))

maxDistance = 0.
bestV = 0.
bestTheta = 0.

for thetaIndex in range(numThetas - 1):
    for velocityIndex in range(numVelocities - 1):
        v0 = velocities[velocityIndex]
        theta0 = thetas[thetaIndex]
        u0 = numpy.array([v0, theta0, x0, y0])# fill 1st element with initial values
        u = u0
        time = 0.
        distance = 0.
        try:
            while time < T and u[3] > 0.:
                u = euler_step(u, f, dt)
                time += dt
            distance = u[2]
        except:
            pass
        # Filter out inf distances
        if distance > maxDistance and distance < 10000:
            maxDistance = distance
            bestV = v0
            bestTheta = theta0
        distances[thetaIndex, velocityIndex] = distance
        
print(maxDistance)
print(bestV)
print(bestTheta)