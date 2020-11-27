import h5py
import matplotlib.pyplot as plt
from matplotlib import cm

### Read cylinder data ###
with h5py.File('dynamic_stall_data.hdf5', 'r') as f:
  rho = f['density'][:] # Density
  rho_u = f['momentum_x'][:] # x-momentum
  rho_v = f['momentum_y'][:] # y-momentum
  p = f['pressure'][:] # Pressure
  x = f['coordinate_x'][:] # x-coordinate
  y = f['coordinate_y'][:] # y-coordinate
  time = f['time'][:] # Time

### Plot a snapshot ###
snap = 120
fig, ax = plt.subplots()
ax.contourf(x, y, rho[:,:,snap], cmap=cm.rainbow)  
ax.set_title('FOM')
ax.axis('equal')
plt.show() 