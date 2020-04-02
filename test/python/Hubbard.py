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
#!/usr/bin/env python

import unittest

from app4triqs import Solver

from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.operators import *
from pytriqs.utility.h5diff import h5diff


class test_hubbard(unittest.TestCase):

    # System Parameters
    U = 1.0
    mu = U / 2
    h = 0.1

    # Construct Parameters
    cp = {}
    cp["beta"] = 10.0
    cp["gf_struct"] = [("up", [0]), ("dn", [0])]
    cp["n_tau"] = 10000
    cp["n_iw"] = 500

    # Set up the Solver
    S = Solver(**cp)
    S.G0_iw["up"] << inverse(iOmega_n + mu + h)
    S.G0_iw["dn"] << inverse(iOmega_n + mu - h)

    # Solve Parameters
    sp = {}
    sp["h_int"] = U * n("up", 0) * n("down", 0)
    sp["max_time"] = -1
    sp["verbosity"] = 3
    sp["post_process"] = True

    # Solve the impurity model
    S.solve(**sp)

    # Store the Result
    with HDFArchive("Hubbard.out.h5", 'w') as arch:
        arch["S"] = S

    # -------- Compare ---------
    # h5diff("hubbard.out.h5", "hubbard.ref.h5")


if __name__ == '__main__':
    unittest.main()
