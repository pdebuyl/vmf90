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

F90 = h5fc -O3 -I.
VPATH = ../src
all : hmf fel
hmf: vmf90_hmf
fel: vmf90_fel

vmf90_version.h: vmf90_version.h.in FORCE
	../scripts/make_version.sh
FORCE:

vmf90.o: vmf90.f90 vmf90_version.h
	$(F90) -c $<

%.o : %.f90
	$(F90) -c $<

vmf90_hmf: h5md.o spline_module.o Vlasov_module.o HMF_module.o ParseText.o vmf90_version.h vmf90.o vmf90_hmf.o
	$(F90) h5md.o spline_module.o Vlasov_module.o HMF_module.o ParseText.o vmf90.o vmf90_hmf.o -o vmf90_hmf

vmf90_fel: h5md.o spline_module.o Vlasov_module.o FEL_module.o ParseText.o vmf90_version.h vmf90.o vmf90_fel.o
	$(F90) h5md.o spline_module.o Vlasov_module.o FEL_module.o ParseText.o vmf90.o vmf90_fel.o -o vmf90_fel

vmf90_ab: h5md.o spline_module.o Vlasov_module.o AB_module.o ParseText.o vmf90_version.h vmf90.o vmf90_ab.o
	$(F90) h5md.o spline_module.o Vlasov_module.o AB_module.o ParseText.o vmf90.o vmf90_ab.o -o vmf90_ab

clean:
	rm -f spline_module.o Vlasov_module.o HMF_module.o AB_module.o ParseText.o vmf90.o FEL_module.o \
	vmf90_hmf.o vmf90_fel.o vmf90_ab.o \
	spline_module.mod vlasov_module.mod parsetext.mod vmf90.mod vmf90.mod \
	hmf_module.mod fel_module.mod ab_module.mod \
	vmf90_version.h vmf90_hmf vmf90_fel \
	h5md.o h5md.mod
