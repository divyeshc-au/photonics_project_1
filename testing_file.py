import numpy as np
import scipy as sp
from scipy.fft import fft2
from scipy.fft import ifft2
from scipy.fft import fftfreq
from scipy.fft import fftshift
import imageio
import cv2

import math
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
    return np.multiply.outer(np.exp( [(((-c*(t)/n + z_)*n*2*np.pi/lam)*1j) for z_ in z ] )*np.heaviside((c*t*1-z),0), J )

def rotation_mat (theta):
    return [[cos(theta),-sin(theta)],[sin(theta),cos(theta)]]

def rotate_mat (matrix, theta):
    return matmul( rotation_mat(-theta), matmul(matrix, rotation_mat(theta)) )

def final_product (n,z,lam,theta,d,J_0):
    j_x = J_0 [0]
    j_y = J_0 [1]
    c = math.cos (theta)
    s = math.sin (theta)

    if z < 0 or z == 0:
        return real_norm (J_0)
    elif  z > 0 and z < d:
        ndarr = (((np.abs (z)+1)>0)*1)
        A = c**2+np.exp((2*np.pi*z/lam)*1j)*(s**2)
        B = c*s*np.exp((2*np.pi*z/lam)*1j)
        C = np.exp((2*np.pi*z/lam)*1j)*(c**2) + s**2
        A = A*ndarr
        B = B*ndarr
        C = C*ndarr
        return (np.real (A*j_x + B*j_y) )**2 + (np.real (B*j_y + C*j_x )**2 )
    else:
        A = c**2+np.exp((2*np.pi*d/lam)*1j)*(s**2)
        B = c*s*np.exp((2*np.pi*d/lam)*1j)
        C = np.exp((2*np.pi*d/lam)*1j)*(c**2) + s**2
        return (np.real (A*j_x + B*j_y) )**2 + (np.real (B*j_y + C*j_x )**2 )

#plt.style.use(['dark_background'])

d = 4
lam = 0.00006 
n = lam/ (2*d)
theta = 0

c = math.cos (theta)
s = math.sin (theta)

x = np.linspace(-10,10,1600) 
xv, yv = np.meshgrid(x, x) 


v = [3,4]
j_0 = np.array(v)/np.linalg.norm([3,4])
j_x = j_0[0]
j_y = j_0[1]

e_field_arr = e_field( j_0,0,xv,lam,n)

#A = c**2+np.exp(complex(0,2*np.pi*d/lam))*(s**2)
#B = c*s*np.exp(complex(0,2*np.pi*d/lam))
#C = np.exp( complex(0,2*np.pi*d/lam))*(c**2) + s**2

#jones_mat_half_wave = np.array [[A,B],[B,C]]

I = final_product (n,xv,lam,theta,d,j_0)
#I = (np.real (A*j_x + B*j_y) )**2 + (np.real (B*j_y + C*j_x )**2 )

#print (e_field( j_0,0,x,lam,n))
#print (e_field (v,0,x,lam,n))
#print([real_norm(E ) for E in  e_field( v,0,x,lam,n)] )
plt.figure(figsize=(5,5))

ax = sb.heatmap (I)

#plt.plot(x, [real_norm( np.matmul( jones_mat_half_wave, E)) for E in  e_field_arr])# ( j_0,0,x,lam,n)] )
#plt.xlabel('Y-Position [mm]')
#plt.ylabel('Z-Position [mm]')
#plt.show()
