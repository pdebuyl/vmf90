#!/usr/bin/env python

# Copyright (C) 2009-2013 Pierre de Buyl

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
  print "        where cmd is plot, snaps or dump and var is expected to be found in"
  print "        the observables of the file. dump will output the time series to"
  print "        standard out."
  exit()

cmd = argv[2]
if (cmd not in ['plot','snaps','dump']):
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

Nx = a['parameters']['Nx'].value
Nv = a['parameters']['Nv'].value
xmin = a['parameters']['xmin'].value
xmax = a['parameters']['xmax'].value
vmin = a['parameters']['vmin'].value
vmax = a['parameters']['vmax'].value
extent = [xmin,xmax,vmin,vmax]

plt.rc('figure.subplot', left = 0.17)
f = plt.figure(figsize=[8.,4.5])

do_show = False

if (cmd == 'plot'):
  n_plot = 0
  for var in argv[3:]:
    if var in a['observables']:
      plt.plot(a['observables'][var]['time'],a['observables'][var]['value'],label=var)
      n_plot+=1
    else:
      print "variable %s not found" % (var,)

  if (n_plot>0):
    plt.xlabel('time')
    plt.legend()
    do_show = True
  else:
    print "available observables: ", a['observables'].keys()

elif (cmd == 'snaps'):
  n = a['fields']['f']['value'].shape[0]
  for i in range(4):
    plt.subplot(2,2,i+1)
    plt.imshow(a['fields']['f']['value'][i*n/4],origin='lower',vmin=0.,extent=extent)
    plt.text(.05,.8,r'$t='+"%.2f" % (a['fields']['f']['time'][i*n/4],)+r'$',transform = plt.gca().transAxes, color='white')
    if i>1: plt.xlabel(r'$x$')
    if i%2==0: plt.ylabel(r'$v$')
  do_show = True

elif (cmd == 'dump'):
  step = None
  var_list = []
  for var in argv[3:]:
    if var in a['observables']:
      if step is None:
        step = a['observables'][var]['step'][:]
        time = a['observables'][var]['time'][:]
      else:
        assert np.allclose(step,a['observables'][var]['step'][:])
      var_list.append(var)
    else:
      print "variable %s not found" % (var,)
      exit()
  if len(var_list) > 0:
    for i in range(len(step)):
      print time[i],
      for var in var_list:
        print "%20.12f " % (a['observables'][var]['value'][i],) ,
      print "\n" ,
  else:
    print "available observables: ", a['observables'].keys()

a.close()

if (do_show): plt.show()
