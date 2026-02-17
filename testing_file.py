import numpy as np
import scipy as sp

import seaborn as sb
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import animation
from matplotlib.animation import PillowWriter


def real_norm (v):
    r_v = np.real (v)
    return r_v[0]**2 + r_v[1]**2

def final_intensity (n,z,lam,theta,d,J_0,t):
    global C
    c = 3* (10**8)/ (1)

    j_x = J_0 [0]*np.exp((((-c*(t)/n + z+10)*n*2*np.pi/lam)*1j))
    j_y = J_0 [1]*np.exp((((-c*(t)/n + z+10)*n*2*np.pi/lam)*1j))
    nd_arr_all1s = (((np.abs (z)+1)>0)*1)
    
    c = math.cos (theta)
    s = math.sin (theta)
    
    A = c**2+np.exp((2*np.pi*(z)*n/lam)*1j)*(s**2)
    B = c*s*np.exp((2*np.pi*(z)*n/lam)*1j)
    C = np.exp((2*np.pi*(z)*n/lam)*1j)*(c**2) + s**2
    
    A2 = c**2+np.exp((2*np.pi*d*n/lam)*1j)*(s**2)
    B2 = c*s*np.exp((2*np.pi*d*n/lam)*1j)
    C2 = np.exp((2*np.pi*d*n/lam)*1j)*(c**2) + s**2
    A2 = A2*nd_arr_all1s
    B2 = B2*nd_arr_all1s
    C2 = C2*nd_arr_all1s
    
    I = nd_arr_all1s*0.0
    I += np.where(z<=0, real_norm (J_0),0)
    I += np.where(np.logical_and((z > 0),(z < d)),  (np.real (A*j_x + B*j_y) )**2 + (np.real (B*j_y + C*j_x )**2 ), 0.0)
    I += np.where(z>=d, (np.real (A2*j_x + B2*j_y) )**2 + (np.real (B2*j_y + C2*j_x )**2), 0.0)
    return I

plt.style.use('dark_background')

lam = 590* (10)**-9 
n = 0.3
d = lam/(2*n)

theta = 0

xlim_min = -1
xlim_max = +1

c = math.cos (theta)
s = math.sin (theta)

x = np.linspace(xlim_min,xlim_max,1600) 
y = np.linspace(xlim_min,xlim_max,1600) 

xv, yv = np.meshgrid(x, y) 


v = [3,4]
j_0 = np.array(v)/np.linalg.norm([3,4])
j_x = j_0[0]
j_y = j_0[1]

I = final_intensity (n,xv,lam,theta,d, j_0,0)

#print (I)

ax = sb.heatmap (I)
ax.set_xlabel ('Z-Position [mm]')
ax.set_ylabel ('Y-Position [mm]')
ax.set_xlim ([xlim_min,xlim_max])
plt.show ()

