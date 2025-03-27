import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def load_trajectory(file):
   # Read the data from the file
   data = np.loadtxt(file, skiprows=1)

   # Separate the time and positions (x, y, z)
   f = {
      't': data[:, 0],        # First column is time
      'pos': data[:, 1:4],    # positions (x, y, z)   
   }
   return f

def load_trajectory_otd(file):
   # Read the data from the file
   data = np.loadtxt(file, skiprows=1)

   # Separate data
   f = {
      't':     data[:, 0],                     # Time
      'bf':    data[:, 1:4],                   # Baseflow (x, y, z)
      'otd1':  data[:, 4:7],                   # OTD basis vector 1 (x, y, z)
      'otd2':  data[:, 7:10],                  # OTD basis vector 2 (x, y, z)
      'cc':    data[:, 10].astype(bool),       # cc flag
      's' :    data[:, 11],                    # sigma
      'sp':    data[:, 12:15],                 # Projected leading eigenvector of Lr_sym
      'EV':    data[:, 15:17].astype(complex), # Instantaneous eigenvalues
      'otd1p': data[:, 17:20],                 # Instantaneous projected OTD mode 1
      'otd2p': data[:, 20:23],                 # Instantaneous projected OTD mode 2
      'FTLE':  data[:, 23:25],                 # Finite-Time Lyapunov Exponents
      'norm':  data[:, 25:27],                 # basis normality
      'orth':  data[:, 27]                     # basis orthogonality
   }
   for i, cc in enumerate(f['cc']):
      if cc:
         pm = [ 1, -1]
         for j in range(2):
            f['EV'][i,j] = f['EV'][i,0] + 1j*pm[j]*f['EV'][i,1]
      else:
         rp = f['EV'][i,:]
         f['EV'][i,:] = sorted(rp, reverse=True)
   if not f['cc'].any():
      f['EV'] = np.real(f['EV'])
   
   return f

def load_LE(file):
    # Read the data from the file
   data = np.loadtxt(file, skiprows=1)

   f = {
      't':      data[:, 0],   # Time
      'period': data[:, 1],   # Period counter
      'LE':     data[:, 2:4], # Lyapunov exponents
      'LE_ref': data[:, 4:6]  # Reference values
   }
   return f

if __name__ == '__main__':
   # Roessler coefficients
   a = 0.2
   b = 0.2
   c = 5.7

   # Extract data
   bf_attr = load_trajectory('roessler_attractor.txt')
   po_otd  = load_trajectory_otd('PO_OTD.txt')
   po_le   = load_LE('PO_LE.txt')
   poc_otd = load_trajectory_otd('PO-chaos_OTD.txt')
   poc_le  = load_LE('PO-chaos_LE.txt')

   # group basis vectors and modes
   basis = (po_otd['otd1'],  po_otd['otd2'])
   modes = (po_otd['otd1p'], po_otd['otd2p'])

   # compute analytical fixed points
   d = np.sqrt(c**2 - 4*a*b)
   fpx = ( ( c - d)/ 2, ( c + d)/ 2 )
   fpy = ( (-c + d)/(2*a), (-c - d)/(2*a) )
   fpz = ( ( c - d)/(2*a), ( c + d)/(2*a) )

   # extract one period
   T = po_le['t'][0]
   npts_p = 0
   while po_otd['t'][npts_p] <= T:
      npts_p += 1
   npts_p -= 1
   # compute cumulative distance along the orbit
   dbf = np.diff(po_otd['bf'][:npts_p], axis=0)
   arclength = np.insert(np.cumsum(np.linalg.norm(dbf, axis=1)), 0, 0)
   arclength /= arclength[-1] # normalize

   # extract chaotic attractor
   rax, ray, raz = bf_attr['pos'].T
   # extract periodic orbit
   pot = po_otd['t']
   pox, poy, poz = po_otd['bf'].T

   #####################################
   #
   #    Plot the roessler attractor 
   #        with periodic orbit
   #
   #####################################

   pox, poy, poz = po_otd['bf'].T
   fig = plt.figure()
   ax = fig.add_subplot(111, projection='3d')
   ax.plot(rax, ray, raz, c='black', alpha=0.5, label='chaotic attractor')
   ax.plot(pox, poy, poz,  c='green', lw=2, label='periodic orbit')
   ax.set_xlabel('x')
   ax.set_ylabel('y')
   ax.set_zlabel('z')
   ax.legend()

   #####################################
   #
   #    Plot the convergence history
   #      of the Lyapunov exponents
   #      computed using OTD modes
   #
   #####################################

   fig, ax = plt.subplots(1, 1, figsize=(10,5))
   for i in range(2):
      ax.plot(po_le['period'], abs(po_le['LE'][:,i]-po_le['LE_ref'][:,i]), 'o-', label=f'Lyapunov exponent {i+1:d}')
   ax.set_yscale('log')
   ax.set_xlabel('periods')
   ax.set_ylabel('Error')
   ax.set_xlim(0, 30)
   ax.legend()
   
   #####################################
   #
   #    Plot the roessler attractor 
   #        with periodic orbit
   #     and the evolution of then
   #          OTD eigenspace
   #
   #####################################

   fig = plt.figure(figsize=(12,6))
   ax = fig.add_subplot(121, projection='3d')
   # Plot the roessler attractor with the periodic orbit
   ax.plot(rax, ray, raz, c='black', alpha=0.5, label='chaotic attractor')
   ax.plot(pox, poy, poz, c='green', lw=2, label='periodic orbit')
   
   # Plot most unstable subspace at 8 points
   npts = 8
   # find equispaced points along orbit
   idx = [ np.argmin( np.abs( arclength - i/npts) ) for i in range(npts) ]
   
   # mark points along the orbit
   ax.scatter(pox[idx], poy[idx], poz[idx], c='black', s=50)
   ax.scatter(pox[0], poy[0], poz[0], s=150, marker='o', edgecolors='k', facecolors='none', label='orbit origin')

   grid_length, vector_l = 3, 3
   # reference x, y meshgrid and 'unit' vectors
   unit = np.linspace(0, grid_length, 2)
   X, Y = np.meshgrid(unit, unit)
   coords = np.stack([X.flatten(),Y.flatten()], axis=1)  # coordinates of the plane [ x | y ]

   # Draw the vectors
   for id in idx:
      # origin of the vectors is the point on the trajectory
      ox, oy, oz = po_otd['bf'][id,:]

      # move the plane to the origin
      Z = np.ones_like(X)*oz
      # prepare surface computation
      base = (X,Y)
      M = np.empty((2,2))
      for i, vec in enumerate(basis):
         M[:,i] = vec[id,:2]
         # add height contribution of the current basis vector v
         Z += base[i]*vec[id,2]

      # compute the surface (linear transformation of the basis vectors (x,y) -> (v1,v2) = M)
      s = np.zeros_like(coords)
      for i, c in enumerate(coords):
         s[i,:] = M @ c
      s += np.outer(np.ones((s.shape[0],)), [ox, oy])
      
      # plot the surface
      ax.plot_surface(s[:,0].reshape(X.shape), s[:,1].reshape(Y.shape), Z, color='black', alpha=0.75, label='OTD subspace')

      # Draw leading singular vector and eigenvectors

      # plot direction of largest instantaneous linear growth
      svx, svy, svz = po_otd['sp'][id,:]
      s1 = np.log(po_otd['s'][id])
      ax.quiver(ox, oy, oz, svx, svy, svz, color='red', length=vector_l, lw=3,  normalize=True, label='OTD largest linear growth')

      for i, mode in enumerate(modes):
         mx, my, mz = mode[id,:]
         ev = np.log(np.abs(np.real(po_otd['EV'][id,i])))
         mplt = ax.quiver(ox, oy, oz, mx, my, mz, color='blue', length=vector_l, lw=2,  normalize=True, label='OTD modes')
   l = 12
   ax.set_xlim(-l, l)
   ax.set_ylim(-l, l)
   ax.set_zlim(0, 2*l)
   ax.set_xlabel('x')
   ax.set_ylabel('y')
   ax.set_zlabel('z')
   handles, labels = plt.gca().get_legend_handles_labels()
   by_label = dict(zip(labels, handles))
   ax.legend(by_label.values(), by_label.keys())

   ax = fig.add_subplot(122)
   # Plot the instantaneous eigenvalues
   plt.axhline(0, color='k')
   plt.axvline(29, color='k', linestyle='--', label='orbit origin')
   for i in range(2):
      ax.plot(po_otd['t']/T, np.real(po_otd['EV'][:,i]), label=r'$\lambda_'+f'{i+1:d}'+'$')
   # Plot sigma max
   ax.plot(po_otd['t']/T, po_otd['s'], label=r'$\sigma_'+f'{{max}}'+'$', color='red', linewidth=2)
   tend = po_otd['t'][-1]
   ax.set_xlabel('period')
   ax.set_ylabel('growth rate')
   ax.set_xlim((tend-2*T)/T,tend/T)
   ax.set_ylim(-4, 4)
   plt.legend()

   plt.show()