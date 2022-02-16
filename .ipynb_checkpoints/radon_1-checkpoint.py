import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import radon_lib as rl


# Testing a time domain forward/adjoint parabolic Radon Transforms
# M. D. Sacchi (March, 2014)
# ----------------------------------------------------------------

Np = 33        # Curvatures
Nt = 250
Nh = 33
dt = 4./1000.
dh = 5. 

p = np.linspace(0.,1.4,Np)
h = np.linspace(0,1000,Nh)

d = np.zeros((Nt,Nh))
m = np.zeros((Nt,Np))

f0 = 34;

wavelet = rl.ricker(dt,f0)

Nw = len(wavelet)
href=h[Nh-1]
m[40:40+Nw,4]=wavelet

m[90:90+Nw,4]=-wavelet

m[75:75+Nw,0]=-wavelet



d = rl.radon_forward(m,Nt,dt,Nh,h,Np,p,href)



m = rl.radon_adjoint(d,Nt,dt,Nh,h,Np,p,href)



mls =rl.radon_cg(d,m*0,Nt,dt,Nh,h,Np,p,href)
 
# Mutes

mls_f =np.copy(mls)
mls_f[:,0:3]=0

dp = rl.radon_forward(mls_f,Nt,dt,Nh,h,Np,p,href)

prim = d-dp

# # Plot results
#==============================================================================

plt.subplot(1,5,1)    
tmp = np.clip(d/np.max(d),-1,1)    
im = plt.imshow(tmp,origin='upper',extent=[min(h),max(h),(Nt-1)*dt,0],cmap=cm.seismic)

plt.axis('tight')
plt.xlabel('Offset (m)')
plt.ylabel('time')
plt.title('Data')


plt.subplot(1,5,2)    
tmp = np.clip(m/np.max(m),-1,1)    
im = plt.imshow(tmp,origin='upper',extent=[min(p),max(p),(Nt-1)*dt,0],cmap=cm.seismic)

plt.axis('tight')
plt.xlabel('Offset (m)')
plt.ylabel('time')
plt.title('Data')

plt.subplot(1,5,3)    
tmp = np.clip(mls/np.max(mls),-1,1)    
im = plt.imshow(tmp,origin='upper',extent=[min(p),max(p),(Nt-1)*dt,0],cmap=cm.seismic)

plt.axis('tight')
plt.xlabel('Offset (m)')
plt.ylabel('time')
plt.title('Data')

plt.subplot(1,5,4)    
tmp = np.clip(prim/np.max(prim),-1,1)    
im = plt.imshow(tmp,origin='upper',extent=[min(h),max(h),(Nt-1)*dt,0],cmap=cm.seismic)


plt.axis('tight')
plt.xlabel('Offset (m)')
plt.ylabel('time')
plt.title('Data')

plt.show()


            