#!/usr/bin/env python

# Copyright (C) 2015 Pierre de Buyl

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

"""Print a configuration file to stdout for a waterbag initial
condition of magnetization M0 and energy U.
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('M0', type=float)
parser.add_argument('U', type=float)
parser.add_argument('--nx', type=int, default=256)
parser.add_argument('--nv', type=int, default=256)
parser.add_argument('--n-moments', type=int, default=0)
parser.add_argument('--n-images', type=int, default=10)
parser.add_argument('--n-top', type=int, default=100)
parser.add_argument('--n-edf', type=int, default=0)
args = parser.parse_args()

import numpy as np
from scipy.optimize import brentq
import os

def magnetization_root(Dtheta, M0):
    """Return zero when Dtheta matches the magnetization."""
    return np.sin(Dtheta)/Dtheta - M0

def solve_for_Dtheta(M0):
    """Return the half-width of the waterbag initial condition Dtheta for
    magnetization M0.
    """
    if M0<1e-12:
        return np.pi
    elif M0>(1.-1e-9):
        return 0.
    else:
        return brentq(magnetization_root, 0.00001, np.pi, args=(M0,))

def HMF_in(values):
    """Writes 'HMF_in' in the directory 'name' with a waterbag initial condition."""
    return """model = HMF
out_file = hmf.h5 
Nx = {nx}
Nv = {nv}
vmin = -3.
vmax = 3.
DT = 0.1
n_steps = 5
n_top = {n_top}
n_images = {n_images}
IC = waterbag
width = {Dth}
bag = {Dp}
Nedf = {n_edf}
Hfield = 0.1
n_moments = {n_moments}
""".format(**values)

Dtheta = solve_for_Dtheta(args.M0)
Dp = np.sqrt(6.*(args.U-0.5*(1.-args.M0**2)))

v = {}
v['Dth'] = Dtheta
v['Dp'] = Dp
v.update(args.__dict__)

print(HMF_in(v))
