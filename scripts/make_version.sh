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

# If by error a "vmf90_version.h" file is in "../src", it will be included
if [ -r ../src/vmf90_version.h ]
then
    echo "Removing src/vmf90_version.h"
    rm ../src/vmf90_version.h
fi

# Check whether we are in a git repository.
if [ "`git rev-parse --is-inside-work-tree 2>/dev/null`" = "true" ]
then
    # Check if the working directory is clean.
    GIT_STATUS="Unclean"
    git status | grep -q "working directory clean" && GIT_STATUS="Clean"
    # Set git_sha1
    GIT_SHA1=`git log -n1 --format="%H"`
    # Set git_describe
    GIT_DESCRIBE=`git describe`
    GIT_DATE=`git log -n1 --format=%aD`
    # Take the original file and set git_status
    sed -e "s/GIT_STATUS/${GIT_STATUS}/g" ../src/vmf90_version.h.in > vmf90_version.h.temp
    sed -i -e "s/GIT_SHA1/${GIT_SHA1}/g" vmf90_version.h.temp
    sed -i -e "s/GIT_DESCRIBE/${GIT_DESCRIBE}/g" vmf90_version.h.temp
    sed -i -e "s/GIT_DATE/${GIT_DATE}/g" vmf90_version.h.temp
else
    # Check for vmf90_version.h.dist
    if [ -r ../src/vmf90_version.h.dist ]
    then
	GIT_DESCRIBE=`cat ../src/vmf90_version.h.dist`
    else
	GIT_DESCRIBE="Out of repository"
    fi
    GIT_STATUS="Out of repository"
    GIT_SHA1="Out of repository"
    GIT_DATE="Out of repository"
    # Take the original file and set git_status
    sed -e "s/GIT_STATUS/${GIT_STATUS}/g" ../src/vmf90_version.h.in > vmf90_version.h.temp
    sed -i -e "s/GIT_SHA1/${GIT_SHA1}/g" vmf90_version.h.temp
    sed -i -e "s/GIT_DESCRIBE/${GIT_DESCRIBE}/g" vmf90_version.h.temp
    sed -i -e "s/GIT_DATE/${GIT_DATE}/g" vmf90_version.h.temp
fi
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
