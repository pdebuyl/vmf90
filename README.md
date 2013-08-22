vmf90 : Vlasov solver for mean-field systems in Fortran 90
==========================================================

Copyright © 2009-2013 Pierre de Buyl

vmf90 is a software for the numerical resolution of the Vlasov equation for 
mean-field systems, currently the Hamiltonian Mean-Field model and the
Colson-Bonifacio model for the free electron laser. vmf90 is based
on the semi-Lagrangian method with cubic spline interpolation.

vmf90 is developed by Pierre de Buyl and is available under the GNU General 
Public License http://www.gnu.org/licenses/gpl.html . The GNU General Public
License is found in the file LICENSE. The homepage for vmf90 is
http://homepages.ulb.ac.be/~pdebuyl/vmf90/ .

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
2. Create a build directory and copy the Makefile:

    mkdir build && cd build
    cp ../scripts/Makefile ./

3. Eventually, edit the ``Makefile`` to adapt to your compiler.
4. Run the make command with the argument ``hmf``:

    make hmf

5. The build directory now contains an executable ``vmf90_hmf``.
6. A ``HMF_in`` configuration file is required to set up the
   simulation. Examples are found in the ``examples`` directory. To copy one of
   these examples, type:

    cp ../scripts/HMF_in.resonances ./HMF_in

7. A run is performed by executing ``vmf90_hmf``:

    ./vmf90_hmf

8. After a run is completed, the data is found in the file ``hmf.h5``. Examples
   scripts on how to read these files in Python are found in the ``scripts``
   directory. To display the total energy and the interaction and kinetic parts,
   issue the following command (NumPy, h5py and Matplotlib are required):

    ../scripts/show_vmf90.py hmf.h5 plot energy en_int kin_kin


Documentation
-------------

A doxygen-generated documentation is readily available by executing:

    doxygen doc/doxy_conf

from the main directory of vmf90. The index of the documentation is found at
``doc/html/index.html``.

Appropriate references to the algorithms are given in the documentation, as well
as a full API documentation for the code and the Fortran modules.

Appropriate citation is appreciated in any scientific publication that makes use
of this software. The latest citation suggestion is:

* P. de Buyl, "vmf90 - A Vlasov solver for mean-field systems in Fortran 90"
    http://homepages.ulb.ac.be/~pdebuyl/vmf90/.
* P. de Buyl, "Numerical resolution of the Vlasov equation for the Hamiltonian
    Mean-Field model", Commun. Nonlinear Sci. Numer. Simulat. vol. 15,
    pp. 2133-2139 (2010).

Please check the vmf90 homepage to see if a more recent reference is appropriate.

