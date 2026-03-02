import matplotlib.pyplot as plt 
import scipy as sp 
import seaborn as sb  
import numpy as np 
from matplotlib.animation import FuncAnimation

global N_X
N_X = 1.49

global N_Y
N_Y = 1.66

global n_x
n_x = 0 

global n_y
n_y = N_Y - N_X 

global lam
lam = 590* (10**-6) # in mm

global alpha
alpha = 1 

global d
d = np.abs (lam/ (2*(N_Y-N_X))) # in mm 

global is_after
is_after = 0

global c
c = 3 * (10*11) # in mm/ s 

global v
v = [1+1j,-1+1j]
v = v/np.linalg.norm (v)

global tot_frame
tot_frame = 200

fig, ax = plt.subplots (figsize= (9,9))
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
#plt.scatter (np.real (v[0]),np.real (v[1]), color='red')
ax.set_xlim (-1,1)
ax.set_ylim (-1,1)

def init ():
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.scatter (np.real (v[0]),np.real (v[1] * np.exp (-1j*np.pi*alpha* is_after)), color='red' )
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
    plt.scatter (np.real (v_updated[0]), np.real (v_updated[1] * np.exp (-1j*np.pi*alpha* is_after) ), color='red')
    ax.set_xlim (-1,1)
    ax.set_ylim (-1,1)


anim = FuncAnimation (fig, update, init_func= init, frames= tot_frame , repeat=False, interval= 33 )
anim.save ('real_vector_before_element.mp4', writer= 'ffmpeg')

is_after = 1
anim = FuncAnimation (fig, update, init_func= init, frames= tot_frame , repeat=False, interval= 33 )
anim.save ('real_vector_after_element.mp4', writer= 'ffmpeg')

alpha = 1.2
anim = FuncAnimation (fig, update, init_func= init, frames= tot_frame , repeat=False, interval= 33 )
anim.save ('real_vector_wrong_wavelength_after_element.mp4', writer= 'ffmpeg')

