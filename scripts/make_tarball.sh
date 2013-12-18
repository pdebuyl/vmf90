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
mv src/vmf90_version.h vmf90/src
sed -i -e '/git_sha1/d' vmf90/src/vmf90_version.h
sed -i -e '/git_status/d' vmf90/src/vmf90_version.h
echo "    character(len=80), parameter :: git_sha1 = 'Out of repository'" >> vmf90/src/vmf90_version.h
echo "    character(len=80), parameter :: git_status = 'Out of repository'" >> vmf90/src/vmf90_version.h
mkdir vmf90/scripts
sed -e '/TOREMOVE/d' scripts/Makefile > vmf90/scripts/Makefile
tar cf TMP.tar vmf90
tar -f ../vmf90-${VMF90_VERSION}.tar --concatenate --update TMP.tar
rm -rf vmf90 TMP.tar
