#!/bin/bash

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

(cd src ; ../scripts/make_version.sh)
VMF90_VERSION=`git describe`
git archive --prefix=vmf90/ -o ../vmf90-${VMF90_VERSION}.tar HEAD
mkdir -p vmf90/src
echo $VMF90_VERSION > vmf90/src/vmf90_version.h.dist
tar -f ../vmf90-${VMF90_VERSION}.tar --update vmf90/
rm -rf vmf90
