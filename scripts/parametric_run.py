#!/usr/bin/env python

# Copyright (C) 2014 Pierre de Buyl

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

"""This program executes a parametric run on the energy U and initial
magnetization M0 for the HMF model. See http://arxiv.org/abs/1310.0805 for
details.
"""

import numpy as np
from scipy.optimize import brentq
import os
import subprocess
from sys import argv, exit

# Parse command line arguments
run=None
if len(argv)>1:
    if argv[1]=='run':
        if len(argv)>2:
            run = True
    elif argv[1]=='show':
        run = False
if run is None:
    print """Usage: %s cmd path_to_vmf90_hmf
       where cmd is run or show. The path to vmf90_hmf is needed fo run only
"""
    exit()

if run:
    assert os.path.isfile(argv[2])
else:
    import matplotlib.pyplot as plt
    import h5py


def magnetization_root(Dtheta, M0):
    """Function returning zero when Dtheta matches the magnetization."""
    return np.sin(Dtheta)/Dtheta - M0

def solve_for_Dtheta(M0):
    """Given a magnetization M0, return the half-width of the waterbag initial
    condition Dtheta
    """
    if M0<1e-12:
        return np.pi
    elif M0>(1.-1e-9):
        return 0.
    else:
        return brentq(magnetization_root, 0.00001, np.pi, args=(M0,))

def write_HMF_in(name, Dth, Dp):
    """Writes 'HMF_in' in the directory 'name' with a waterbag initial condition."""
    f = file(os.path.join(name,'HMF_in'), 'w')
    f.write("""model = HMFext
out_file = hmf.h5 
Nx = 256
Nv = 256
vmin = -3.
vmax = 3.
DT = 0.1
n_steps = 10
n_top = 40
n_images = 5
IC = waterbag
width = %f 
bag = %f 
Nedf = 0 
Hfield = 0.1
""" % (Dth,Dp))
    f.close()

U_data = []
M0_data = []
av_data = []
for U in np.linspace(0.51,0.71,11):
    for M0 in np.linspace(0.,1., 11):
        name = "PR_%05.3f_%05.3f" % (U, M0)
        Dtheta = solve_for_Dtheta(M0)
        Dp = np.sqrt(6.*(U-0.5*(1.-M0**2)))
        if run:
            os.mkdir(name)
            write_HMF_in(name, Dtheta, Dp)
            subprocess.call([argv[2]], cwd=name)
        else:
            f = h5py.File(os.path.join(name,'hmf.h5'), 'r')
            mask = f['observables/Mx/time'].value >= 20
            U_data.append(U)
            M0_data.append(M0)
            av_data.append(f['observables/Mx/value'][mask].mean())

if not run:
    plt.figure()
    plt.scatter(M0_data, U_data, s=800, c=av_data)
    plt.xlim(0, 1)
    plt.ylim(0.5, 0.71)
    plt.xlabel(r'$M_0$')
    plt.ylabel(r'$U$')
    plt.title(r'Average of $m_x$')
    plt.colorbar()
    plt.show()

