#!/usr/bin/env python

# Copyright (C) 2009-2011 Pierre de Buyl

# This file is part of vmf90

# vmf90 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# vmf90 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with vmf90.  If not, see <http://www.gnu.org/licenses/>.

from sys import argv

if (len(argv)<3):
  print "Usage : %s show_vmf90.py filename.h5 cmd var1 [var2] [...]" % argv[0]
  print "        where cmd is plot or snaps and var is expected to be found in the"
  print "        observables of the file."
  exit()

cmd = argv[2]
if (cmd not in ['plot','snaps']):
  print "Command %s unknown" % (cmd,)
  exit()

import h5py

try:
  a = h5py.File(argv[1],'r')
except:
  print "Problems opening file", argv[1]
  exit()

import numpy as np
import matplotlib.pyplot as plt

d = dict()
i=0
for n in a['info/time_names'].value:
  d[n] = i
  i+=1

names = a['data'].keys()
n = len(names)
t = a['time'].value
Nx = a['info/Nx'].value
Nv = a['info/Nv'].value
xmin = a['info/xmin'].value
xmax = a['info/xmax'].value
vmin = a['info/vmin'].value
vmax = a['info/vmax'].value
extent = [xmin,xmax,vmin,vmax]

plt.rc('figure.subplot', left = 0.17)
f = plt.figure(figsize=[10.,5.2])

if (cmd == 'plot'):
  n_plot = 0
  for var in argv[3:]:
    if (var in d):
      idx = d[var]
      plt.plot(t[0,:],t[idx,:],label=var)
      n_plot+=1
    else:
      print "variable %s not found" % (var,)

  if (n_plot>0):
    plt.legend()
  else:
    print "available observables: ", a['info/time_names'].value

elif (cmd == 'snaps'):
  for i in range(3):
    plt.subplot(1,3,i+1)
    plt.imshow(a['data'][names[i*n/3]]['f'],origin='lower',vmin=0.,extent=extent)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$v$')


a.close()

plt.show()

