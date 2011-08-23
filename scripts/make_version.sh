#!/usr/bin/bash
# Check if the working directory is clean.
GIT_STATUS=`git status --porcelain`
if [ -z "$GIT_STATUS" ]
then
GIT_STATUS="Clean"
else
GIT_STATUS="Unclean"
fi
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
