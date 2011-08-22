F90 = h5fc -O3
VPATH = ../src
all : runHMF

%.o : %.f90
	$(F90) -c $<

runHMF: spline_module.o Vlasov_module.o HMF_module.o ParseText.o runHMF.o
	$(F90) spline_module.o Vlasov_module.o HMF_module.o ParseText.o runHMF.o -o runHMF

clean:
	rm -f spline_module.o Vlasov_module.o HMF_module.o ParseText.o runHMF.o \
	spline_module.mod Vlasov_module.mod HMF_module.mod ParseText.mod runHMF.mod \
	runHMF