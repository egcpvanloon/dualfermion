/*******************************************************************************
 *
 * app4triqs: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
 *
 * app4triqs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * app4triqs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * app4triqs. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include <app4triqs/solver_core.hpp>

#include <triqs/gfs.hpp>
#include <triqs/h5.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace app4triqs;

TEST(app4triqs, HubbardAtom) { // NOLINT

  // System Parameters
  double U  = 1.0;
  double mu = U / 2;
  double h  = 0.1;

  // Construct Parameters
  constr_params_t cp;
  cp.beta      = 10.0;
  cp.gf_struct = {{"up", {0}}, {"dn", {0}}};
  cp.n_tau     = 10000;
  cp.n_iw      = 500;

  // Set up the Solver
  solver_core S(cp);
  int up = 0, dn = 1;
  S.G0_iw[up](iw_) << 1.0 / (iw_ + mu + h);
  S.G0_iw[dn](iw_) << 1.0 / (iw_ + mu - h);

  // Solve Parameters
  solve_params_t sp;
  sp.h_int        = U * n("up", 0) * n("down", 0);
  sp.max_time     = -1;
  sp.verbosity    = 3;
  sp.post_process = true;

  // Solve the impurity model
  S.solve(sp);

  // Store the Result
  {
    auto arch = triqs::h5::file("hubbard.out.h5", 'w');
    h5_write(arch, "S", S);
  }

  // Compare against the reference data
  // h5diff("hubbard.out.h5", "hubbard.ref.h5")
}

int main(int argc, char **argv) {                                                                                                                  \
  ::mpi::environment env(argc, argv);                                                                                                              \
  ::testing::InitGoogleTest(&argc, argv);                                                                                                          \
  return RUN_ALL_TESTS();                                                                                                                          \
}
