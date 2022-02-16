# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 20:31:49 2015

@author: msacchi
"""
import numpy as np
import matplotlib.pyplot as plt

def wigb(d,dt,h,xcur,color):

# Plot wiggle seismic plot (python version of Xin-gong Li faumous wigb.m)

    [nt,nx] = np.shape(d)
    dmax = np.max(d)
    d = d/dmax
    t = np.linspace(0,(nt-1)*dt,nt)
    tmax = np.amax(t)
    hmin = np.amin(h)
    hmax = np.amax(h)
   
    c = xcur*np.mean(np.diff(h))

    plt.axis([hmin-2*c, hmax+2*c, tmax, 0.])
    d[nt-1,:]=0
    d[0,:]=0
    for k in range(0,nx):
        s =d[:,k]*c
        plt.plot(s+h[k], t, color,linewidth=1)
        b = h[k]+s.clip(min=0) 
        plt.fill(b,t,color)

    return
    
def ricker(dt,f0):    
         
#Ricker wavelet of central frequency f0 sampled every dt seconds

# M D Sacchi, 2015,  Email: msacchi@ualberta.ca

    nw = 2.5/f0/dt
    nw = 2*int(nw/2)
    nc = int(nw/2)
    a = f0*dt*3.14159265359
    n = a*np.arange(-nc,nc)
    b = n**2
    return (1-2*b)*np.exp(-b)


    


