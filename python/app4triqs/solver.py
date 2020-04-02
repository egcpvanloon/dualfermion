###############################################################################
#
# app4triqs: A TRIQS based impurity solver
#
# Copyright (c) 2019 The Simons foundation
#   authors: Nils Wentzell
#
# app4triqs is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# app4triqs is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# app4triqs. If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################
from solver_core import SolverCore

from pytriqs.gf import *
from pytriqs.utility import mpi


# === The SolverCore Wrapper

class Solver(SolverCore):
    def __init__(self, beta, gf_struct, n_iw=500, n_tau=5001):
        """
        Initialise the solver.

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
        n_iw : integer [default=500]
               Number of Matsubara frequencies used for the Green's functions.
        n_tau : integer [default=5001]
               Number of imaginary time points used for the Green's functions.
        """
        # Initialise the core solver
        SolverCore.__init__(self,
                            beta=beta,
                            gf_struct=gf_struct,
                            n_iw=n_iw,
                            n_tau=n_tau)
        self.gf_struct = gf_struct
        self.n_iw = n_iw
        self.n_tau = n_tau

    def solve(self, **params_kw):
        """
        Solve the impurity problem.

        Parameters
        ----------
        params_kw : dict {'param':value} that is passed to the core solver.
                     The only two required parameters are
                        * `h_int`: The local interaction Hamiltonian
                        * `n_cycles`: The number of Monte-Carlo cycles
                     For the other optional parameters see documentation.
                     Note that in this Python Wrapper the alpha-tensor is optional.
                     If not given, it will be constructed from the density matrix of
                     the SC Hartree Fock solution.
        """

        h_int = params_kw['h_int']
        gf_struct = self.gf_struct

        # Call the core solver's solve routine
        return SolverCore.solve(self, **params_kw)
