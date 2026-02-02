import numpy as np
import scipy as sp
from scipy.fft import fft2
from scipy.fft import ifft2
from scipy.fft import fftfreq
from scipy.fft import fftshift
import imageio
import cv2

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import animation
from matplotlib.animation import PillowWriter
import pint
import scienceplots

u = pint.UnitRegistry() 
c = 3* (10**8)/ (1)

def real_norm (v):
    return np.linalg.norm( np.real( v))

def e_field (J,t,z,lam,n):
    phi = ((-c*(t)/n + z)*n*2*np.pi/lam)
    return np.multiply.outer(np.exp( [complex( 0,((-c*(t)/n + z_)*n*2*np.pi/lam)) for z_ in z ] )*np.heaviside(-(c*t*1-z),0), J )

def rotation_mat (theta):
    return [[cos(theta),-sin(theta)],[sin(theta),cos(theta)]]

def rotate_mat (matrix, theta):
    return matmul( rotation_mat(-theta), matmul(matrix, rotation_mat(theta)) )

#plt.style.use(['dark_background'])

d = 4
lam = 0.00000660 
n = lam/ (2*d)

x = np.linspace(-2,2,1600) 
xv, yv = np.meshgrid(x, x) 

U0 = (np.abs(xv)< d/2) * (np.abs(yv)<0.5) # This will create the double slit 
U0 = U0.astype(float) 

v = [3,4]
j_0 = np.array(v)/np.linalg.norm([3,4])

jones_mat_half_wave = np.array([[1,0],[0,np.exp( complex(0,2*np.pi*d/lam))]])

#print (e_field( j_0,0,x,lam,n))
#print (e_field (v,0,x,lam,n))
#print([real_norm(E ) for E in  e_field( v,0,x,lam,n)] )
plt.figure(figsize=(5,5))
plt.plot(x, [real_norm( np.matmul( jones_mat_half_wave, E)) for E in  e_field( j_0,0,x,lam,n)] )
plt.xlabel('Y-Position [mm]')
plt.ylabel('Z-Position [mm]')
plt.show()
