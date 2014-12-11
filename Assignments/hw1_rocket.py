# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 07:52:29 2014

@author: Gary
"""

import numpy
#import math

m_s = 50.
g = 9.81
rho = 1.091
A = 3.1415927 * 0.5 * 0.5
v_e = 325.
C_D = 0.15
m_p0 = 100.
deltaT = 0.1
burnRate = 20.

def f(u, time):
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

    h = u[0]    
    v = u[1]
    m_p = max(m_p0 - (burnRate * time), 0)
    m_dot_p = burnRate
    if time >= 4.99:
        m_dot_p = 0
    
    m = m_s + m_p
    drag = 0.5 * rho * v * abs(v) * A * C_D
   
    return numpy.array([v,
                        (1./m) * ((-m * g) + (m_dot_p * v_e) - drag),
                        m_p,
                        m_dot_p])


def euler_step(u, f, dt, time):
    return u + dt * f(u, time)

T = 50.
t = 0.
N = int(T / deltaT)
times = numpy.linspace(0., T, N)

# initialize the array containing the solution
u = numpy.empty((N, 4))
u[0] = numpy.array([0., 0., m_p0, burnRate])

foundMaxHeight = False
foundZero = False

maxVelocity = 0.
maxVelocityTime = 0.
maxVelocityHeight = 0.

# time loop
for n in range(N-1):
    newU = euler_step(u[n], f, deltaT, t)         ### call euler_step() ###
    t += deltaT
    
    u[n+1] = newU
    
    if newU[1] <= 0 and not foundMaxHeight:
        foundMaxHeight = True
        print "foundMaxHeight at time %f prevHeight:%f prevVel:%f height:%f vel:%f" % (t, u[n][0], u[n][1], u[n+1][0], u[n+1][1])
        
    if newU[0] <= 0 and t > deltaT and not foundZero:
        foundZero = True
        print "returned to zero at time %f prevHeight:%f prevVel:%f height:%f vel:%f" % (t, u[n][0], u[n][1], u[n+1][0], u[n+1][1])
        
    if newU[1] > maxVelocity:
        maxVelocity = newU[1]
        maxVelocityTime = t
        maxVelocityHeight = newU[0]
        
velocities = u[:,1]
maxvelocity = max(velocities)
print "Max Velocity at %f val %f %f height:%f" % (maxVelocityTime, maxVelocity, maxvelocity, maxVelocityHeight)

