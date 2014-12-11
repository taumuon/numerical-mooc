# -*- coding: utf-8 -*-
"""
Created on Thu Oct 02 19:23:57 2014

@author: Gary
"""

import numpy
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# have a solution for rho against x and t
# Velocity can be obtained for any given rho by
def calcVFromRho(rho):
    return V_max * (1. - (rho / rho_max))
    
def calcFFromRho(rho):
    return rho * calcVFromRho(rho)

def stepRho(nt, dt, nx, rhoInit):
    return numpy.ones(nx)
    
def fPrime(rho, rho_max, V_max):
    return V_max - (2. * V_max * rho / rho_max)

def convertKmhToMs(kmh):
    return kmh * 1000 / 60. / 60.

#V_max = 80.
V_max = 136.
L = 11.
rho_max = 250.
nx = 51
dt = 0.001
dx = L/(nx - 1)

x = numpy.linspace(0, L, nx)
rho0 = numpy.ones(nx) * 20
rho0[10:20] = 50

plt.plot(x, rho0)

V = numpy.ones(nx)
for i in range(nx):
    V[i] = calcVFromRho(rho0[i])

# plt.plot(x, V)

minV = numpy.min(V)
minVConverted = convertKmhToMs(minV)
print minVConverted

# boundary condition:
# rho (0, t) = 10

rho = rho0.copy()

# want three minutes. dt = 0.001 hours
t_desired = 0.05 # 3 minutes in fraction of an hour
nt = 50 # t_desired / dt

rho_boundary = 20.

for n in range(1,nt):
  rhon = rho.copy() 
  for i in range(1,nx): 
    rho[1:] = rhon[1:]-fPrime(rhon[1:], rho_max, V_max)*dt/dx*(rhon[1:]-rhon[0:-1])
    rho[0] = rho_boundary

V_threeMin = numpy.ones(nx)
for i in range(nx):
    V_threeMin[i] = calcVFromRho(rho[i])

V_threeMinAverage = numpy.average(V_threeMin)
print convertKmhToMs(V_threeMinAverage)

V_threeMinMin = numpy.min(V_threeMin)
print convertKmhToMs(V_threeMinMin)