F90 = h5fc -O3
VPATH = ../src
all : runHMF

vmf90_version.h: vmf90_version.h.in FORCE
	../scripts/make_version.sh
FORCE:

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
