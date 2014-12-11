import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

uvinitial = numpy.load('./uvinitial.npz')
U = uvinitial['U']
V = uvinitial['V']

#fig = plt.figure(figsize=(8,5))
#plt.subplot(121)
#plt.imshow(U, cmap = cm.RdBu)
#plt.xticks([]), plt.yticks([]);
#plt.subplot(122)
#plt.imshow(V, cmap = cm.RdBu)
#plt.xticks([]), plt.yticks([]);

n = 192
Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.065 # Bacteria 1

dh = 5./(n-1)

#T = 8000
T = 8000

dt = .9 * dh**2 / (4.*max(Du, Dv))

nt = int(T/dt)

def ftcs(U, V, nt, Du, Dv, dt, dx, dy):
    for n in range(nt):
        Un = U.copy()
        Vn = V.copy()
        U[1:-1,1:-1] = Un[1:-1,1:-1] + Du *\
            (dt/dy**2 * (Un[2:,1:-1] - 2*Un[1:-1,1:-1] + Un[:-2,1:-1]) +\
             dt/dx**2 * (Un[1:-1,2:] - 2*Un[1:-1,1:-1] + Un[1:-1,:-2])) -\
             dt * (Un[1:-1,1:-1]*Vn[1:-1,1:-1]*Vn[1:-1,1:-1]) +\
             dt * (F-(F * Un[1:-1,1:-1]))
  
        V[1:-1,1:-1] = Vn[1:-1,1:-1] + Dv *\
            (dt/dy**2 * (Vn[2:,1:-1] - 2*Vn[1:-1,1:-1] + Vn[:-2,1:-1]) +\
             dt/dx**2 * (Vn[1:-1,2:] - 2*Vn[1:-1,1:-1] + Vn[1:-1,:-2])) +\
             dt * (Un[1:-1,1:-1]*Vn[1:-1,1:-1]*Vn[1:-1,1:-1]) -\
             dt *(F+k)*Vn[1:-1,1:-1]
  
        # Enforce Neumann BCs
        U[:,-1] = U[:,-2]
        U[-1,:] = U[-2,:]
        U[:,0] = U[:,1]
        U[0,:] = U[1,:]

        V[:,-1] = V[:,-2]        
        V[-1,:] = V[-2,:]
        V[:,0] = V[:,1]
        V[0,:] = V[1,:]
    return U

Unew = ftcs(U, V, nt, Du, Dv, dt, dh, dh)

fig = plt.figure(figsize=(8,5))
plt.subplot(121)
plt.imshow(Unew, cmap = cm.RdBu)
plt.xticks([]), plt.yticks([]);
plt.subplot(122)
plt.imshow(Vnew, cmap = cm.RdBu)
plt.xticks([]), plt.yticks([]);

ans = Unew[100,::40]
print ans