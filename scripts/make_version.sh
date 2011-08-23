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

# Check if the working directory is clean.
GIT_STATUS="Unclean"
git status | grep -q "working directory clean" && GIT_STATUS="Clean"
# Take the original file and set git_status
sed -e "s/GIT_STATUS/${GIT_STATUS}/g" ../src/vmf90_version.h.in > vmf90_version.h.temp
# Set git_sha1
GIT_SHA1=`git log -n1 --format="%H"`
sed -i -e "s/GIT_SHA1/${GIT_SHA1}/g" vmf90_version.h.temp
# Set git_describe
GIT_DESCRIBE=`git describe`
GIT_DATE=`git log -n1 --format=%aD`
GIT_DESCRIBE="${GIT_DESCRIBE} \/ ${GIT_DATE}"
sed -i -e "s/GIT_DESCRIBE/${GIT_DESCRIBE}/g" vmf90_version.h.temp
# If there is already a file vmf90_version.h, replace it only if the newly
# generated file is different, to prevent the triggering of the makefile rule
if [ -r vmf90_version.h ]
then
! cmp vmf90_version.h.temp vmf90_version.h && mv vmf90_version.h.temp vmf90_version.h
else
mv vmf90_version.h.temp vmf90_version.h
fi
# Remove the temporary file if it is still there
rm -f vmf90_version.h.temp
