import numpy
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

nx = 81
dx = .25
dt = 0.0002
gamma = 1.4
L = 20.

t = 0.01
# nt = 50 #t / dt
nt = 50

rho_l = 1.
rho_r = 0.125
u_l = 0.
u_r = 0.
p_l = 100000.
p_r = 10000.

rho_0 = numpy.ones(nx) * rho_l
rho_0[40:81] = rho_r

u_0 = numpy.ones(nx) * u_l
u_0[40:81] = u_r

p_0 = numpy.ones(nx) * p_l
p_0[40:81] = p_r

x = numpy.linspace(0, L, nx)

#plt.plot(x, p_0)

uvec_0 = numpy.empty((nx, 3))
for initIndex in range(0, nx):
    e = p_0[initIndex] / ((gamma - 1.) * rho_0[initIndex])
    e_T = e + (0.5 * u_0[initIndex] * u_0[initIndex])
    uvec_0[initIndex][0] = rho_0[initIndex]
    uvec_0[initIndex][1] = rho_0[initIndex] * u_0[initIndex]
    uvec_0[initIndex][2] = (rho_0[initIndex] * e_T)


def get_f_from_u(uvec_i):
    fnew = numpy.ones(3)
    a = (uvec_i[1] * uvec_i[1] / uvec_i[0])
    fnew[0] = uvec_i[1]
    fnew[1] = a + ((gamma - 1.) * (uvec_i[2] - (0.5 * a)))
    fnew[2] = (uvec_i[1] / uvec_i[0]) * (uvec_i[2] + ((gamma - 1.) * (uvec_i[2] - 0.5*a)))
    return fnew

uvec = numpy.empty((nx, 3))
for ii in range(0, nx):
    uvec[ii][0] = uvec_0[ii][0]
    uvec[ii][1] = uvec_0[ii][1]
    uvec[ii][2] = uvec_0[ii][2]

for n in range(0, nt):
    uvec_n = numpy.empty((nx, 3))
    for ii in range(0, nx):
        uvec_n[ii][:] = uvec[ii][:]

    fvec_n = numpy.empty((nx, 3))
    for ii in range(0, nx):
        fvec_n[ii] = get_f_from_u(uvec_n[ii])

    #Richtmyer
    uvec_n_plus_half = numpy.empty((nx - 1, 3))
    for i in range(0, nx - 1):
        uvec_n_plus_half[i] =  0.5 * (uvec_n[i + 1][:] + uvec_n[i][:]) - ((dt / (2. * dx)) * (fvec_n[i + 1][:] - fvec_n[i][:]))

    fvec_n_plus_half = numpy.empty((nx - 1, 3))
    for i in range(0, nx - 1):
        fvec_n_plus_half[i] = get_f_from_u(uvec_n_plus_half[i])
    for i in range(1, nx - 1):
        #uvec_i_n_plus_one = numpy.ones(3)
        #uvec_i_n_plus_one[:] = 
        uvec[i][:] = uvec_n[i][:] - ((dt / dx)* (fvec_n_plus_half[i] - fvec_n_plus_half[i - 1]))
        
finalRho = numpy.ones(nx)
for ii in range (0, nx):
    finalRho[ii] = uvec[ii][0]
    
finalU = numpy.ones(nx)
for ii in range (0, nx):
    finalU[ii] = uvec[ii][1]/uvec[ii][0]

finalP = numpy.ones(nx)
for ii in range (0, nx):
    finalP[ii] = (gamma - 1) * (uvec[ii][2] - (0.5*(uvec[ii][1]*uvec[ii][1]/uvec[ii][0])))

xmap = numpy.linspace(0, L, nx)
for xi in range (0, nx):
    xmap[xi] = x[xi] - 10.
    
plt.plot(xmap, finalU)
# plt.plot(x, finalRho)

index = 50
velocityAt2_5 = finalU[index]
pressureAt2_5 = finalP[index]
densityAt2_5 = finalRho[index]
print velocityAt2_5
print pressureAt2_5
print densityAt2_5
print xmap[40]
print xmap[50]