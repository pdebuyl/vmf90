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

F90 = h5fc -O3
VPATH = ../src
all : runHMF

vmf90_version.h: vmf90_version.h.in FORCE    # TOREMOVE
	../scripts/make_version.sh           # TOREMOVE
FORCE:                                       # TOREMOVE

vmf90.o: vmf90.f90 vmf90_version.h
	$(F90) -c $<

%.o : %.f90
	$(F90) -c $<

runHMF: spline_module.o Vlasov_module.o HMF_module.o ParseText.o vmf90_version.h vmf90.o runHMF.o
	$(F90) spline_module.o Vlasov_module.o HMF_module.o ParseText.o vmf90.o runHMF.o -o runHMF

clean:
	rm -f spline_module.o Vlasov_module.o HMF_module.o ParseText.o vmf90.o runHMF.o \
	spline_module.mod Vlasov_module.mod HMF_module.mod ParseText.mod vmf90.mod runHMF.mod \
	vmf90_version.h runHMF
