import matplotlib.pyplot as plt 
import scipy as sp 
import seaborn as sb  
import numpy as np 

plt.style.use ('dark_background')

N_X = 1.49
N_Y = 1.66

n_x = 0 
n_y = N_Y - N_X 

lam = 590* (10**-6) # in mm
d = np.abs (lam/ (2*(N_Y-N_X))) # in mm 

global c
c = 3 * (10*11) # in mm/ s 

lim_left = -(10**-3)
lim_right = (10**-2)

z = np.linspace (lim_left,lim_right,1000)
y = np.linspace (lim_left,lim_right,1000)

Z,Y = np.meshgrid (z,y)

def I_1 (z,v,lamd):
    v_x = v [0]* np.exp (1j*2*np.pi*z/lamd)
    v_y = v [1]* np.exp (1j*2*np.pi*z/lamd)
    return np.real (v_x)**2 + np.real (v_y)**2
def I_2 (z,v,phi_x,phi_y):
    v_x = v [0]
    v_y = v [1]
    return np.real (v_x*np.exp (phi_x*z))**2 + np.real (v_y*np.exp (phi_y*z))**2 
def I_3 (z,v,d,lamd,phi_x,phi_y):
    v_x = v [0]
    v_y = v [1]
    return np.real (v_x* np.exp ( 1j*2*np.pi*(z-d)/lamd + phi_x*d ))**2 + np.real (v_y* np.exp ( 1j*2*np.pi*(z-d)/lamd + phi_y*d ))**2


def intensity (z,n_x,n_y,lamd,d,v):
    v_x = v [0]
    v_y = v [1]
    phi_x = 1j*2*np.pi*n_x/lamd
    phi_y = 1j*2*np.pi*n_y/lamd
    return (z<=0)*I_1 (z,v,lamd) + np.logical_and (z>0,z<d)*I_2 (z,v,phi_x,phi_y) + (z>=d)* I_3 (z,v,d,lamd,phi_x,phi_y) 


def intensity_plot (alpha=1,beta=0,v=[3,4]):
    intsty = intensity ( -c*beta+Z, n_x, n_y, lam *alpha , d , v)
    ax = sb.heatmap (intsty)
    ax.invert_yaxis ()
    plt.show ()

def update_plot (val):
    a = alpha.val
    b = beta.val 
    intensity_plot (a,b,[3,4])

v = [3,4]

intensity_plot ()

