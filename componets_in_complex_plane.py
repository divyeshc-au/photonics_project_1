import matplotlib.pyplot as plt 
import scipy as sp 
import seaborn as sb  
import numpy as np 
from matplotlib.animation import FuncAnimation

N_X = 1.49
N_Y = 1.66

n_x = 0 
n_y = N_Y - N_X 

lam = 590* (10**-6) # in mm
d = np.abs (lam/ (2*(N_Y-N_X))) # in mm 

global c
c = 3 * (10*11) # in mm/ s 

global v
v = np.array([1+1j,-1+1j])
v = v/np.linalg.norm (v)

global alpha 
alpha = 1 

global is_after 
is_after = 0

global tot_frame
tot_frame = 200

fig, ax = plt.subplots (figsize= (9,9))
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
#plt.scatter (np.real (v[0]),np.imag (v[0]), color='red', label= 'E_x')
#plt.scatter (np.real (v[1]),np.imag (v[1]), color='blue', label= 'E_y')
#plt.legend ()
ax.set_xlim (-1,1)
ax.set_ylim (-1,1)
def init ():
    v_init = v
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.scatter (np.real (v_init[0]),np.imag (v_init[0]), color='red', label= 'E_x')
    plt.scatter (np.real (v_init[1]*np.exp (-1j*np.pi*alpha*is_after)
),np.imag (v_init[1]*np.exp (-1j*np.pi*alpha*is_after)
), color='blue', label= 'E_y')
    plt.legend ()
    ax.set_xlim (-1,1)
    ax.set_ylim (-1,1)

def update (frame):
    v_updated = v* np.exp (1j*2*np.pi*frame/tot_frame)
    ax.clear ()
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.scatter (np.real (v_updated[0]),np.imag (v_updated[0]), color='red', label= 'E_x')
    plt.scatter (np.real (v_updated[1]*np.exp (-1j*np.pi*alpha*is_after)
),np.imag (v_updated[1]*np.exp (-1j*np.pi*alpha*is_after)
), color='blue', label= 'E_y')
    plt.legend ()
    ax.set_xlim (-1,1)
    ax.set_ylim (-1,1)


anim = FuncAnimation (fig, update, init_func= init, frames= tot_frame , repeat=False, interval= 33 )
anim.save ('vector_in_complex_plane_before_element.mp4', writer= 'ffmpeg')

is_after = 1 
anim = FuncAnimation (fig, update, init_func= init, frames= tot_frame , repeat=False, interval= 33 )
anim.save ('vector_in_complex_plane_after_element.mp4', writer= 'ffmpeg')

alpha = 1.2
anim = FuncAnimation (fig, update, init_func= init, frames= tot_frame , repeat=False, interval= 33 )
anim.save ('vector_in_complex_plane_wrong_wavelength.mp4', writer= 'ffmpeg')

