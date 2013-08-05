Configuration files syntax {#config_file_syntax}
==========================

The vmf90 programs vmf90_hmf and vmf90_fel read the parameters for a simulation
in the files HMF_in and FEL_in, respectively.

Example parameter files are found in the ``scripts`` directory.

General syntax
--------------

Parameters are given according to the syntax:

    variable_name = variable_value


Parameters common to all vmf90 are requested to setup the simulation:

``Nx``
    Number of grid points in x-direction.

``Nv``
    Number of grid poits in v-direction.

``vmax``
    Upper boundary of the box in the v-direction.

``vmin``
    Lower boundary of the box in the v-direction. Not need for HMF_in where vmin
    is taken as -vmax.

``DT``
    Time step for the inner loop. This is the elementary timestep for the
    algorithm.

``n_steps``
    Number of repetitions of the inner loop. The observables are not computed
    within the inner loop.

``n_top``
    Number of repetitions of the outer loop. At each iteration, the observables
    are computed and written to the output file. The total number of timesteps
    is n_steps*n_top and the total time of the simulation is n_steps*n_top*DT.

``n_images``
    Number snapshots of phase space that are written to disk.

``IC``
    Initial condition. This is program-specific.


vmf90_hmf syntax
----------------

vmf90_hmf supports the following extra options:

``Nedf``
    Number of data points for storing the energy distribution function.

``Hfield``
    For fixed-field simulations, the value of the field.

Initial conditions for vmf90_hmf
--------------------------------

A waterbag. Calls \ref vlasov_module::init_carre

    IC = waterbag
    width = angular half width
    bag = velocity half width

A waterbag with a cosine perturbation. Calls \ref vlasov_module::init_carre

    IC = wb_eps
    width = angular half width
    bag = velocity half width
    epsilon = small perturbation
