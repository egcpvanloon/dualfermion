################################################################################
#
# dualfermion: second-order dual fermion perturbation theory
#
# Copyright (C) 2020, E.G.C.P. van Loon
# 
# Based on TRIQS: a Toolbox for Research in Interacting Quantum Systems
# 
# dualfermion is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# dualfermion is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# dualfermion. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

from dpt_core  import DptCore
from pytriqs.gf import *
import pytriqs.utility.mpi as mpi
import numpy as np

class Dpt(DptCore):
#TODO: which parameters to give at initialization, which at run?

    def __init__(self, beta, gf_struct, n_iw=128, n_iw2=32,n_iW=31,N_x=1,N_y=1,N_z=1):
        """
        Initialise the dual perturbation theory.

        Parameters
        ----------
        beta : scalar
               Inverse temperature.
        gf_struct : list of pairs [ [str,[int,...]], ...]
                    Structure of the Green's functions. It must be a
                    list of pairs, each containing the name of the
                    Green's function block as a string and a list of integer
                    indices.
                    For example: ``[ ['up', [0, 1, 2]], ['down', [0, 1, 2]] ]``.
        n_iw : integer, optional
               Number of Matsubara frequencies used for the Green's functions.
        n_iw2 : integer, optional
               Number of Fermionic Matsubara frequencies used for the Two-Particle Green's functions.
        n_iw : integer, optional
               Number of Bosonic Matsubara frequencies used for the Two-Particle Green's functions.
        """
        if isinstance(gf_struct,dict):
            print "WARNING: gf_struct should be a list of pairs [ [str,[int,...]], ...], not a dict"
            gf_struct = [ [k, v] for k, v in gf_struct.iteritems() ]

        # Initialise the dual perturbation theory. 
        DptCore.__init__(self, beta=beta, gf_struct=gf_struct, 
                            n_iw=n_iw,n_iw2=n_iw2,n_iW=n_iW,N_x=N_x,N_y=N_y,N_z=N_z)

        self.gf_struct = gf_struct
        self.n_iw = n_iw
        self.n_iw2= n_iw2
        self.n_iW = n_iW
        self.N_x  = N_x
        self.N_y  = N_y
        self.N_z  = N_z

    def run(self, **params_kw):
        """
        Run the dual perturbation theory once.
        Large parts of the interface with the core program run via hdf5 files

        Parameters
        ----------
        params_kw : dict {'param':value} that is passed to the core program.
        """

        # Pre-processing
        params_kw['sigmad_subset'] = params_kw.get('sigmad_subset',self.gf_struct)

        # Call the core solver's solve routine
        status = DptCore.run(self, **params_kw)

        # Post-processing:

        return status
