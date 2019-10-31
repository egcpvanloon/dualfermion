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
#include "./solver_core.hpp"
#include "./post_process.hpp"

namespace app4triqs {

  solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

    // Initialize the non-interacting Green function
    G0_iw = block_gf<imfreq>{{p.beta, Fermion, p.n_iw}, p.gf_struct};

    // Initialize the result containers
    G_tau    = block_gf<imtime>{{p.beta, Fermion, p.n_tau}, p.gf_struct};
    G_iw     = G0_iw;
    Sigma_iw = G0_iw;
  }

  // -------------------------------------------------------------------------------

  void solver_core::solve(solve_params_t const &solve_params) {

    last_solve_params = solve_params;

    if (world.rank() == 0)
      std::cout << "\n"
                   "APP4TRIQS Solver\n";

    // Assert hermiticity of the given Weiss field
    if (!is_gf_hermitian(G0_iw)) TRIQS_RUNTIME_ERROR << "Please make sure that G0_iw fullfills the hermiticity relation G_ij[iw] = G_ji[-iw]*";

    // Merge constr_params and solve_params
    params_t params(constr_params, solve_params);

    // Reset the results
    container_set::operator=(container_set{});

    // TODO Solve the impurity model

    // Post Processing
    if (params.post_process) { post_process(params); }
  }

  // -------------------------------------------------------------------------------

  void solver_core::post_process(params_t const &p) {

    if (world.rank() == 0)
      std::cout << "\n"
                   "Post-processing ... \n";

    // TODO
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write_attribute(grp, "TRIQS_HDF5_data_scheme", solver_core::hdf5_scheme());
    h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(AS_STRING(TRIQS_GIT_HASH)));
    h5_write_attribute(grp, "APP4TRIQS_GIT_HASH", std::string(AS_STRING(APP4TRIQS_GIT_HASH)));
    h5_write(grp, "", s.result_set());
    h5_write(grp, "constr_params", s.constr_params);
    h5_write(grp, "last_solve_params", s.last_solve_params);
    h5_write(grp, "G0_iw", s.G0_iw);
  }

  solver_core solver_core::h5_read_construct(triqs::h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "", s.result_set());
    h5_read(grp, "last_solve_params", s.last_solve_params);
    h5_read(grp, "G0_iw", s.G0_iw);
    return s;
  }

} // namespace app4triqs
