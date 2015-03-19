vmf90 : Vlasov solver for mean-field systems in Fortran 90
==========================================================

Copyright © 2009-2015 Pierre de Buyl

vmf90 is a software for the numerical resolution of the Vlasov equation for 
mean-field systems, currently the Hamiltonian Mean-Field model and the
Colson-Bonifacio model for the free electron laser. vmf90 is based
on the semi-Lagrangian method with cubic spline interpolation.

vmf90 is developed by Pierre de Buyl and is available under the [GNU General 
Public License](http://www.gnu.org/licenses/gpl.html).
The GNU General Public License is found in the file LICENSE.
The homepage for vmf90 is <https://github.com/pdebuyl/vmf90>.

vmf90 is presented in P. de Buyl, [The vmf90 program for the numerical
resolution of the Vlasov equation for mean-field systems
](http://dx.doi.org/10.1016/j.cpc.2014.03.004), Comp. Phys. Comm. (2014) -
[[arXiv.org:1310.0805]](http://arxiv.org/abs/1310.0805).
Citations to this reference are recommended and appreciated if you use vmf90 to
obtain scientific results.


Requirements
------------

* A Fortran 95 compiler.
* The [HDF5](http://www.hdfgroup.org/HDF5/) library.

To analyze the data, Python example scripts are provided. They require

* [Python](http://python.org/).
* The [NumPy](http://numpy.scipy.org/) library.
* The [Matplotlib](http://matplotlib.sourceforge.net/) library.
* The [h5py](http://www.h5py.org/) library.

Usage
-----

1. Make sure that the ``h5fc`` Fortran compiler script for HDF5 is properly set
   up.
2. Run the make command with the argument ``hmf``:

    make hmf

3. The build directory now contains an executable ``vmf90_hmf``.
4. A ``HMF_in`` configuration file is required to set up the
   simulation. Examples are found in the ``scripts`` directory. To copy one of
   these examples, type:

    cd build
    cp ../scripts/HMF_in.resonances ./HMF_in

7. A run is performed by executing ``vmf90_hmf``:

    ./vmf90_hmf

8. After a run is completed, the data is found in the file ``hmf.h5``. Examples
   scripts on how to read these files in Python are found in the ``scripts``
   directory. To display the total energy and the interaction and kinetic parts,
   issue the following command (NumPy, h5py and Matplotlib are required):

    ../scripts/show_vmf90.py hmf.h5 plot energy en_int en_kin

To use the FEL program, the instructions are similar but "make hmf" becomes
"make fel". Also, the program for the FEL is vmf90_fel and the parameters file
for the FEL program is FEL_in.

Documentation
-------------

A doxygen-generated documentation is readily available by executing:

    make doc

from the main directory of vmf90. The index of the documentation is found at
``doc/html/index.html``.

Appropriate references to the algorithms are given in the documentation, as well
as a full API documentation for the code and the Fortran modules.

