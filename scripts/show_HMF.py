#!/usr/bin/env python
from sys import argv

if (len(argv)<3):
  print "Usage : %s show_HMF.py filename.h5 cmd" % argv[0]
  print "        where cmd is one of en, mass, m, mxy, rho, phi, snaps"
  exit()

import h5py
import numpy as np
import matplotlib.pyplot as plt

a = h5py.File(argv[1],'r')

t = a['time'].value
names = a['data'].keys()
n = len(names)
Nx = a['info/Nx'].value
Nv = a['info/Nv'].value
xmin = a['info/xmin'].value
xmax = a['info/xmax'].value
vmin = a['info/vmin'].value
vmax = a['info/vmax'].value
extent = [xmin,xmax,vmin,vmax]

plt.rc('figure.subplot', left = 0.17)
f = plt.figure(figsize=[10.,5.2])

if (argv[2]=='en'):
  for i in [2,3,4]:
    plt.plot(t[0,:],t[i,:])
  plt.legend(['en','int','kin'])
  plt.xlabel(r'$t$',fontsize=18)
elif (argv[2]=='mass'):
  plt.plot(t[0,:],t[1,:])
  plt.legend([r'$mass$'])
  plt.xlabel(r'$t$',fontsize=20)
elif (argv[2]=='m'):
  plt.plot(t[0,:],np.sqrt(t[6,:]**2+t[7,:]**2))
  plt.legend([r'$m$'])
  plt.xlabel(r'$t$',fontsize=18)
elif (argv[2]=='mxy'):
  plt.plot(t[0,:],t[6,:])
  plt.plot(t[0,:],t[7,:])
  plt.legend([r'$m_x$',r'$m_y$'])
  plt.xlabel(r'$t$',fontsize=18)
elif (argv[2]=='snaps'):
  for i in range(3):
    plt.subplot(1,3,i+1)
    plt.imshow(a['data'][names[i*n/3]]['f'],origin='lower',vmin=0.)
elif (argv[2]=='rho'):
  plt.plot(np.linspace(xmin,xmax,Nx),a['data'][names[-1]]['rho'])
  plt.xlabel(r'$\theta$',fontsize=18)
  plt.ylabel(r'$\rho(\theta)$',fontsize=18)
elif (argv[2]=='phi'):
  plt.plot(np.linspace(vmin,vmax,Nv),a['data'][names[-1]]['phi'])
  plt.xlabel(r'$p$',fontsize=18)
  plt.ylabel(r'$\phi(p)$',fontsize=18)

a.close()

plt.show()

