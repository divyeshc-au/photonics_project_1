import matplotlib.pyplot as plt 
import scipy as sp 
import seaborn as sb  
import numpy as np 
from matplotlib.animation import FuncAnimation

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

v = [3,4]

def I_1 (z,v,lamd,t):
    v_x = v [0]* np.exp (1j*2*np.pi*(c*t-z)/lamd)
    v_y = v [1]* np.exp (1j*2*np.pi*(c*t-z)/lamd)
    return np.real (v_x)**2 + np.real (v_y)**2

def I_2 (z,v,lamd,phi_x,phi_y,t):
    v_x = v [0]* np.exp (1j*2*np.pi*c*t/lamd-phi_x*z)
    v_y = v [1]* np.exp (1j*2*np.pi*c*t/lamd-phi_y*z)
    return np.real (v_x)**2 + np.real (v_y)**2 

def I_3 (z,v,d,lamd,phi_x,phi_y,t):
    v_x = v [0]* np.exp ( 1j*2*np.pi*(c*t-z+d)/lamd - phi_x*d )
    v_y = v [1]* np.exp ( 1j*2*np.pi*(c*t-z+d)/lamd - phi_y*d )
    return np.real (v_x)**2 + np.real (v_y)**2

def intensity (z,n_x,n_y,lamd,d,v,t):
    phi_x = 1j*2*np.pi*n_x/lamd
    phi_y = 1j*2*np.pi*n_y/lamd
    return (z<=0)*I_1 (z,v,lamd,t) + np.logical_and (z>0,z<d)*I_2 (z,v,lamd,phi_x,phi_y,t) + (z>=d)* I_3 (z,v,d,lamd,phi_x,phi_y,t) 

tot_frame = 50
fig, ax = plt.subplots (figsize= (16,9))
sb.heatmap (intensity_matrices [0])

def init_heatmap ():
    I = intensity ( Z, n_x, n_y, lam , d , v, 0)
    sb.heatmap ( I, cbar=False )

def update_heatmap (frame):
    I = intensity ( Z, n_x, n_y, lam , d , v, frame) 
    ax.clear ()
    sb.heatmap ( I, cbar=False)

#def intensity_plot (alpha=1,v=[3,4],t=0):
#    intsty = intensity ( Z, n_x, n_y, lam *alpha , d , v, t)
#    ax = sb.heatmap (intsty)
#    ax.invert_yaxis ()
#    plt.show ()

#for i in range (5):
#    intensity_plot (t= i)

#intensity_plot ()

anim = FuncAnimation (fig, update_heatmap, init_func= init_heatmap, frames= tot_frame , repeat=False, interval= 16 )
#plt.show ()
anim.save ('heatmap_test_animation.mp4', writer= 'ffmpeg')

