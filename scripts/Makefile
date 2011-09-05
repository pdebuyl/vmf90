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
all : hmf fel
hmf: vmf90_hmf
fel: vmf90_fel

vmf90_version.h: vmf90_version.h.in FORCE    # TOREMOVE
	../scripts/make_version.sh           # TOREMOVE
FORCE:                                       # TOREMOVE

vmf90.o: vmf90.f90 vmf90_version.h
	$(F90) -c $<

%.o : %.f90
	$(F90) -c $<

vmf90_hmf: spline_module.o Vlasov_module.o HMF_module.o ParseText.o vmf90_version.h vmf90.o vmf90_hmf.o
	$(F90) spline_module.o Vlasov_module.o HMF_module.o ParseText.o vmf90.o vmf90_hmf.o -o vmf90_hmf

vmf90_fel: spline_module.o Vlasov_module.o FEL_module.o ParseText.o vmf90_version.h vmf90.o vmf90_fel.o
	$(F90) spline_module.o Vlasov_module.o FEL_module.o ParseText.o vmf90.o vmf90_fel.o -o vmf90_fel

clean:
	rm -f spline_module.o Vlasov_module.o HMF_module.o ParseText.o vmf90.o vmf90_hmf.o vmf90_fel.o \
	spline_module.mod vlasov_module.mod hmf_module.mod parsetext.mod vmf90.mod \
	FEL_module.o fel_module.mod \
	vmf90_version.h vmf90_hmf vmf90_fel
