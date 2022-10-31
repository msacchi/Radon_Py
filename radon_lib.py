
"""
Created on Mon Feb 16 12:00:56 2015

@author: msacchi
"""

from numba import jit
import numpy as np

@jit(nopython=True)
def radon_adjoint(d,Nt,dt,Nh,h,Np,p,href):

# Adjoint Time-domain Parabolic Radon Operator

# d(nt,nh): data
# dt      : sampling interval
# h(Nh)   : offset
# p(Np)   : curvature of parabola
# href    : reference offset
# Returns m(nt,np) the Radon coefficients 

# M D Sacchi, 2015,  Email: msacchi@ualberta.ca

 
    m=np.zeros((Nt,Np))

    for itau in range(0,Nt):
        for ih in range(0,Nh):
            for ip in range(0,Np):
                t = (itau)*dt+p[ip]*(h[ih]/href)**2
                it = np.int(t/dt)
                if it<Nt:
                    if it>0:
                        m[itau,ip] +=  d[it,ih]
    
    return m
        

def radon_forward(m,Nt,dt,Nh,h,Np,p,href):

# Forward Time-domain Parabolic Radon Transform

# m(nt,nh): Radon coefficients 
# dt      : sampling interval
# h(Nh)   : offset
# p(Np)   : curvature of parabola
# href    : reference offset
# Returns d(nt,nh) the synthetized data from the Radon coefficients

# M D Sacchi, 2015,  Email: msacchi@ualberta.ca

    d=np.zeros((Nt,Nh))
      
    for itau in range(0,Nt):
        for ih in range(0,Nh):
            for ip in range(0,Np):
                t = (itau)*dt+p[ip]*(h[ih]/href)**2
                it=np.int(t/dt)
                if it<Nt:
                    if it>=0:
                        d[it,ih] +=  m[itau,ip]                   
    return d
    

def radon_cg(d,m0,Nt,dt,Nh,h,Np,p,href,Niter):
         
# LS Radon transform. Finds the Radon coefficients by minimizing
# ||L m - d||_2^2 where L is the forward Parabolic Radon Operator.
# The solution is found via CGLS with operators L and L^T applied on the
# flight

# M D Sacchi, 2015,  Email: msacchi@ualberta.ca

    m = m0  
    
    s = d-radon_forward(m,Nt,dt,Nh,h,Np,p,href) # d - Lm
    pp = radon_adjoint(s,Nt,dt,Nh,h,Np,p,href)  # pp = L's 
    r = pp
    q = radon_forward(pp,Nt,dt,Nh,h,Np,p,href)
    old = np.sum(np.sum(r*r))
    print("iter","  res")

    for k in range(0,Niter):
         alpha = np.sum(np.sum(r*r))/np.sum(np.sum(q*q))
         m +=  alpha*pp
         s -=  alpha*q
         r = radon_adjoint(s,Nt,dt,Nh,h,Np,p,href)  # r= L's
         new = np.sum(np.sum(r*r))
         print(k, new)
         beta = new/old
         old = new
         pp = r + beta*pp
         q = radon_forward(pp,Nt,dt,Nh,h,Np,p,href) # q=L pp
           
    return m 


